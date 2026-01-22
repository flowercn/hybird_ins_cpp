#ifndef INS_WRAPPER_H
#define INS_WRAPPER_H

#include "ins_state.h"
#include "kf_state.h"
#include "glv.h"
#include "earth.h"
#include "ins_math.h" 
#include <vector>
#include <iostream>

struct AlignResult {
    Eigen::Vector3d att; 
    Eigen::Vector3d vn;  
    Eigen::Vector3d pos; 
    Eigen::Vector3d eb;  
    Eigen::Vector3d db;  
    bool success;
};

struct IMUData { 
    Eigen::Vector3d wm, vm; 
    double t; 
};

struct AlignConfig {
    Eigen::Vector3d att_ref;      
    Eigen::Vector3d pos_ref;      
    Eigen::Vector3d phi_init_err; 
    Eigen::Vector3d wvn_err;      
    double eb_sigma;   
    double db_sigma;   
    double web_psd;    
    double wdb_psd;    
};

// [修改] 将原函数重构为增量式对准引擎
class AlignmentEngine {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // 内部状态
    KFAlignVN kf;
    Earth eth;
    Eigen::Quaterniond qnb;
    Eigen::Vector3d vn;
    Eigen::Vector3d pos;
    
    // 缓冲逻辑 (用于 nn=2)
    int step_counter;
    int nn;
    double nts;
    IMUData last_data; // 缓存上一帧
    bool finished;

    AlignmentEngine() : finished(false) {}

    // 初始化 (对应原函数的前半部分)
    void Init(double ts, const GLV& glv, const AlignConfig& config) {
        eth = Earth(glv);
        eth.eupdate(config.pos_ref, Eigen::Vector3d::Zero()); 
        
        vn = Eigen::Vector3d::Zero();
        pos = config.pos_ref; 
        qnb = INSMath::a2qua(config.att_ref);
        
        // 设定两子样逻辑
        nn = 2; 
        nts = nn * ts;
        step_counter = 0;
        
        // KF 初始化
        kf.Init(nts, glv, 
                0.0,                   
                config.web_psd,        
                config.wdb_psd,        
                config.eb_sigma,       
                config.db_sigma,       
                config.wvn_err(0));    

        kf.Pxk.block<3,3>(0,0) = config.phi_init_err.array().square().matrix().asDiagonal();
        finished = false;
        
        std::cout << "[AlignEngine] Initialized. Step size nn=" << nn << std::endl;
    }

    // [核心] 单步执行：每次 global_T 增加 ts 时调用一次
    // 只有当攒够 nn 帧时，内部才会真正执行 KF 更新
    void Step(const IMUData& curr_data) {
        using namespace Eigen;

        // 1. 如果是两子样的第一帧，仅缓存
        if (step_counter == 0) {
            last_data = curr_data;
            step_counter++;
            return; 
        }

        // 2. 如果是两子样的第二帧，开始计算
        if (step_counter == 1) {
            const auto& d1 = last_data;
            const auto& d2 = curr_data; // current
            
            // --- 严格复刻原逻辑 ---
            Vector3d wmm = d1.wm + d2.wm;
            Vector3d vmm = d1.vm + d2.vm;
            
            // 1. 圆锥误差 (Coning)
            Vector3d dphim = 2.0/3.0 * d1.wm.cross(d2.wm);
            Vector3d phim = wmm + dphim; 
            
            // 2. 划船误差 (Sculling)
            Vector3d scullm = 2.0/3.0 * (d1.wm.cross(d2.vm) + d1.vm.cross(d2.wm));
            
            // 3. 旋转误差
            Vector3d rotm = 0.5 * wmm.cross(vmm);
            
            Vector3d dvbm = vmm + rotm + scullm;

            // --- 机械编排 ---
            Matrix3d Cnn = INSMath::rv2m(-eth.wien * nts / 2.0); // 注意这里的 Cnn 近似
            Matrix3d Cnb = INSMath::q2mat(qnb); 
            Vector3d dvn = Cnn * Cnb * dvbm;
            vn = vn + dvn + eth.gn * nts;
            
            Vector3d rv_phi = phim;
            Vector3d rv_zeta = eth.winn * nts; 
            qnb = INSMath::qupdt2(qnb, rv_phi, rv_zeta);

            // --- KF 更新 ---
            kf.UpdatePhi(dvn, Cnb, eth.wien, nts);
            kf.Update(vn); // 这里假设 vn 为零速观测，或者外部传入的 vn_ref
            
            // 反馈校正
            kf.Feedback(qnb, vn);

            // 重置计数器
            step_counter = 0;
        }
    }

    // 获取结果用于切换到 INS
    AlignResult GetResult() {
        AlignResult res;
        res.att = INSMath::m2att(INSMath::q2mat(qnb)); 
        res.vn  = vn;
        res.pos = pos; 
        res.eb = kf.xk.segment<3>(6); 
        res.db = kf.xk.segment<3>(9); 
        res.success = true;
        return res;
    }
};

#endif // INS_WRAPPER_H