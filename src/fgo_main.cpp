#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <Eigen/Dense>

#include "support.h"

using namespace std;
using namespace Eigen;

int main() {
    // // 1. 基础配置
    // GLV glv; 
    // string fog_path = "../fog.csv"; 
    // double ts = 1.0 / 400.0;
    
    // // 加载数据 (support.h 内部已包含自动计算 scale_ratio 的逻辑)
    // auto all_data = load_and_process_data(fog_path, ts, glv.g0);
    // if (all_data.empty()) return -1;

    // // 2. 配置参数 (与 MATLAB 代码对齐)
    // AlignConfig align_cfg;
    // align_cfg.att_ref << 0.03404385, 0.3160598, 0.6178784; // 这里的数值根据你的MATLAB代码调整
    // align_cfg.att_ref *= glv.deg;  
    
    // align_cfg.pos_ref << 32.0286 * glv.deg, 118.8533 * glv.deg, 17.0; 
    
    // align_cfg.phi_init_err << 0.01, 0.01, 0.1; align_cfg.phi_init_err *= glv.deg; 
    // align_cfg.wvn_err << 0.001, 0.001, 0.001;
    
    // // 器件误差设置
    // align_cfg.eb_sigma = 0.2 * glv.dph;          
    // align_cfg.db_sigma = 100 * glv.ug;           
    // align_cfg.web_psd  = 0.021 * glv.dpsh;       
    // align_cfg.wdb_psd  = 10 * glv.ugpsHz;

    // // ---------------------------------------------------------
    // // 步骤 A: 运行标准对准 (获取 Cnb 和 Bias)
    // // ---------------------------------------------------------
    // double t_align = 300.0;
    // SinsEngine align_engine(ts);
    // align_engine.Init(align_cfg, 30.0, 270.0); // Coarse 30s + Fine 270s = 300s

    // cout << "------------------------------------------------" << endl;
    // cout << "Step 1: Running Initial Alignment (300s)..." << endl;
    
    // for (const auto& d : all_data) {
    //     if (d.t > t_align) break; // 只跑前300秒
    //     align_engine.Step(d.wm, d.vm, d.t);
    // }

    // // 获取对准结果
    // Matrix3d Cnb_aligned = INSMath::q2mat(align_engine.GetQnb());
    // Vector3d eb_est = align_engine.GetBiasGyro();
    // Vector3d db_est = align_engine.GetBiasAcc();
    
    // cout << "[Align Result] T=300s" << endl;
    // cout << "  Att (deg): " << (INSMath::m2att(Cnb_aligned) * glv.rad).transpose() << endl;
    // cout << "  EB (deg/h):" << (eb_est * glv.rad * 3600.0).transpose() << endl;
    // cout << "  DB (ug):   " << (db_est / glv.ug).transpose() << endl;

    // // ---------------------------------------------------------
    // // 步骤 B: 数据回补 (Rotation Compensation)
    // // ---------------------------------------------------------
    // cout << "\nStep 2: Rotating raw data to N-Frame (0,0,0)..." << endl;
    
    // vector<IMUData> rotated_data;
    // rotated_data.reserve(all_data.size());

    // for (const auto& d : all_data) {
    //     IMUData d_new = d;
    //     // [关键逻辑] 将 b 系增量左乘 Cnb，投影到 n 系
    //     // 这样新的 d_new 看起来就像是一个水平放置(姿态为0)的 IMU 输出的数据
    //     d_new.wm = Cnb_aligned * d.wm; 
    //     d_new.vm = Cnb_aligned * d.vm;
    //     rotated_data.push_back(d_new);
    // }

    // // ---------------------------------------------------------
    // // 步骤 C: 零姿态纯惯导验证
    // // ---------------------------------------------------------
    // cout << "\nStep 3: Running Zero-Attitude Pure INS Verification..." << endl;

    // // 1. 初始化纯导状态
    // // 注意：这里姿态强制设为 0 (Identity Matrix)
    // Vector3d att_zero = Vector3d::Zero(); 
    // Vector3d vn_zero  = Vector3d::Zero();
    // Earth eth(glv);
    // eth.eupdate(align_cfg.pos_ref, vn_zero);

    // INSState ins_verify(att_zero, vn_zero, align_cfg.pos_ref, ts, eth);
    
    // // 2. 注入刚才对准估出来的零偏 (补偿掉，看残差)
    // ins_verify.set_bias(eb_est, db_est);
    // ins_verify.vertical_damping_mode = 1; // 开启高程阻尼，方便看水平发散

    // // 3. 准备记录
    // string out_name = "res_rotated_ins.csv";
    // ofstream out_file(out_name);
    // out_file << "t,lat,lon,h,vn,ve,vd,roll,pitch,yaw,bg_x,bg_y,bg_z,ba_x,ba_y,ba_z,status" << endl;
    // out_file << fixed << setprecision(9);

    // int log_cnt = 0;
    // for (const auto& d : rotated_data) {
        
    //     // 纯惯导更新
    //     ins_verify.update(d.wm, d.vm, glv, eth);

    //     // 记录
    //     if (++log_cnt % 40 == 0) { // 10Hz 记录
    //         Vector3d att_deg = ins_verify.att * glv.rad;
    //         Vector3d pos_deg = ins_verify.pos; 
    //         pos_deg(0) *= glv.rad; pos_deg(1) *= glv.rad;
            
    //         // 此时的 eb/db 是固定的补偿值
    //         Vector3d eb_show = ins_verify.eb * glv.rad * 3600.0;
    //         Vector3d db_show = ins_verify.db / glv.ug;

    //         out_file << d.t << ","
    //                  << pos_deg(0) << "," << pos_deg(1) << "," << pos_deg(2) << ","
    //                  << ins_verify.vn(0) << "," << ins_verify.vn(1) << "," << ins_verify.vn(2) << ","
    //                  << att_deg(1) << "," << att_deg(0) << "," << att_deg(2) << "," // R,P,Y
    //                  << eb_show(0) << "," << eb_show(1) << "," << eb_show(2) << ","
    //                  << db_show(0) << "," << db_show(1) << "," << db_show(2) << ","
    //                  << 3 << endl; // Status=3 (Nav)
    //     }
    // }

    // cout << "Finished. Results saved to " << out_name << endl;
    // return 0;
}