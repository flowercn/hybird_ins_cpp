#ifndef SINS_ENGINE_H
#define SINS_ENGINE_H

#include <vector>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "glv.h"
#include "earth.h"
#include "ins_math.h"
#include "ins_state.h"
#include "kf_state.h"

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

enum class EngineState {
    IDLE,            
    COARSE_ALIGNING, 
    FINE_ALIGNING,   
    NAVIGATING,      
    FAULT            
};

class SinsEngine {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    GLV glv;
    Earth eth;
    INSState ins;       
    EngineState state;

    AlignConfig align_cfg;  
    double coarse_duration; 
    double fine_duration;   

    // Coarse 变量
    Eigen::Vector3d coarse_acc_wm; 
    Eigen::Vector3d coarse_acc_vm; 
    int coarse_sample_cnt;         

    // Fine/KF 变量
    KFAlignVN align_kf;            
    Eigen::Quaterniond align_qnb;  
    Eigen::Vector3d align_vn;      
    Eigen::Vector3d align_pos;     
    
    // 缓冲逻辑
    int align_step_counter;
    int align_nn;
    double align_nts;
    IMUData align_last_data; 
    bool align_finished;

    // --- [新增] 迭代对准核心变量 ---
    int kf_round;                  // 当前是第几轮 KF (1 或 2)
    double time_kf_switch;         // 第一轮结束、切换第二轮的时间点
    double scale_ratio;            // 自动计算的加速度标度系数
    Eigen::Vector3d accum_bias_gyro; // 第一轮估计并固化的陀螺零偏
    Eigen::Vector3d accum_bias_acc;  // 第一轮估计并固化的加计零偏

    SinsEngine(double ts = 1.0/400.0);
    void Init(const AlignConfig& cfg, double coarse_time_s, double fine_time_s);
    void Step(const Eigen::Vector3d& wm, const Eigen::Vector3d& vm, double t);

    // Getters
    Eigen::Vector3d GetAttDeg() const; 
    Eigen::Vector3d GetVel() const;    
    Eigen::Vector3d GetPosDeg() const; 
    Eigen::Vector3d GetBiasGyro() const; 
    Eigen::Vector3d GetBiasAcc() const;  
    Eigen::Quaterniond GetQnb() const;
    void InjectBias(const Eigen::Vector3d& db_gyro, const Eigen::Vector3d& db_acc);

    static Eigen::Vector3d AlignCoarse(const Eigen::Vector3d& wmm, const Eigen::Vector3d& vmm, double latitude);
};

#endif // SINS_ENGINE_H