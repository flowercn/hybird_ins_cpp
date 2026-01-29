/**
 * @file test_atom_feedback.cpp
 * @brief 原子陀螺辅助 FOG 静态导航测试 - 实时反馈模式
 * * 逻辑：
 * 1. 执行混合对准，锁定初始姿态 Cnb_ref。
 * 2. 计算参考基准 wie_b_ref = Cnb_ref' * wie_n。
 * 3. 进入导航循环：
 * - 每 2秒 (T_cycle) 累积 FOG 数据。
 * - 模拟原子陀螺测量 (基于 wie_b_ref + 原子噪声)。
 * - 计算 FOG 零偏误差: bias_est = fog_mean - atom_meas。
 * - 将估计的零偏注入 INS，修正后续解算。
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "sins_engine.h"
#include "cai_sim.h"
#include "support.h"

using namespace std;
using namespace Eigen;
using namespace cai;

int main() {
    // 1. 加载数据
    double ts = 1.0 / 400.0;
    GLV glv;
    Vector3d pos_ref(32.0286 * glv.deg, 118.8533 * glv.deg, 17.0);
    Earth eth(glv); eth.update(pos_ref, Vector3d::Zero());
    double local_g = eth.gn.norm(); // 获取当地重力模长

    auto imu_data = LoadIMUData("../fog3h.csv", ts, local_g);
    cout << "Loaded " << imu_data.size() << " epochs." << endl;

    SinsEngine ins(ts);
    ins.Sins_Init(pos_ref, Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero());
    
    HybridAlignConfig cfg;
    cfg.t_coarse = 60;
    cfg.t_fine = 3600; // 1小时精对准
    cfg.verbose = false;
    
    auto align_res = ins.Run_HybridAlign(imu_data, cfg);
    if (!align_res.valid) { cerr << "Alignment failed!" << endl; return -1; }
    
    cout << "Alignment done. Initial Att: " << align_res.att.transpose() /glv.deg << " deg" << endl;

    // 3. 初始化原子陀螺仿真器
    AtomicGyroSimulator atom(pos_ref, glv);
    atom.Init(align_res.att); // 锁定 reference wie_b
    atom.PrintConfig();

    // 4. 准备导航 (从对准结束时刻开始)
    size_t nav_start_idx = static_cast<size_t>(cfg.t_fine / ts);
    
    ins.eth.update(pos_ref, Vector3d::Zero());
    ins.ins = INSState(align_res.att, Vector3d::Zero(), pos_ref, ts, ins.eth);
    // 这里的 db 是 0，因为 Run_HybridAlign 还原后没有估计加计零偏
    ins.ins.set_bias(align_res.eb, align_res.db); 

    cout << "Initial Bias (deg/h): " << (align_res.eb / glv.deg * 3600.0).transpose() << endl;
    cout << "Initial Bias (ug): " << (align_res.db * 1e6).transpose() << endl;
    
    //exit(0);
    // 窗口参数
    double T_cycle = 2.0;
    int samples_per_cycle = static_cast<int>(T_cycle / ts);
    Vector3d fog_sum = Vector3d::Zero();
    int sample_count = 0;
    
    // [New] 加速度计统计变量
    Vector3d acc_sum_nav = Vector3d::Zero();
    long nav_samples = 0;

    // 记录日志
    ofstream log("atom_feedback_nav.csv");
    log << "time,lat_err_m,lon_err_m,h_err_m,drift_m,eb_x,eb_y,eb_z\n";

    cout << "\n[3] Navigation with Atomic Feedback..." << endl;
    double max_drift = 0;
    
    for (size_t i = nav_start_idx; i < imu_data.size(); ++i) {
        // [New] 累积导航过程中的加速度计读数
        acc_sum_nav += imu_data[i].vm / ts; // 累积比力 (m/s^2)
        nav_samples++;
        
        // 4.1 累积 FOG 数据
        fog_sum += imu_data[i].wm / ts; // 转换为角速率 rad/s
        sample_count++;

        // 4.2 执行单步导航
        ins.Step_Nav(imu_data[i].wm, imu_data[i].vm);

        // 4.3 周期性反馈 (原子陀螺)
        if (sample_count >= samples_per_cycle) {
            //差分估计
            Vector3d fog_mean = fog_sum / sample_count;
            Vector3d atom_meas = atom.Measure(); 
            Vector3d bias_new = fog_mean - atom_meas;
            ins.ins.set_bias(bias_new, align_res.db);
            // 重置计数器
            fog_sum.setZero();
            sample_count = 0;
        }

        // 4.4 记录与统计
        if ((i - nav_start_idx) % 400 == 0) { // 1Hz log
            double lat_err = (ins.ins.pos(0) - pos_ref(0)) * glv.Re;
            double lon_err = (ins.ins.pos(1) - pos_ref(1)) * glv.Re * cos(pos_ref(0));
            double h_err = ins.ins.pos(2) - pos_ref(2);
            double drift = sqrt(lat_err*lat_err + lon_err*lon_err);
            if(drift > max_drift) max_drift = drift;

            Vector3d eb_dph = ins.ins.eb / glv.deg * 3600.0;
            
            log << fixed << setprecision(4) 
                << (i - nav_start_idx) * ts << "," 
                << lat_err << "," << lon_err << "," << h_err << "," << drift << ","
                << eb_dph(0) << "," << eb_dph(1) << "," << eb_dph(2) << "\n";
        }
    }
    
    log.close();
    
    // 5. 结果输出
    double total_time_h = (imu_data.size() - nav_start_idx) * ts / 3600.0;
    cout << "--------------------------------" << endl;
    cout << "Total Nav Time: " << total_time_h << " h" << endl;
    cout << "Max Drift:      " << max_drift << " m" << endl;
    cout << "Final Drift:    " << max_drift / 1852.0 << " nm" << endl;
    
    // [New] 加速度计均值分析
    cout << "\n[Diagnostic] Accelerometer Statistics:" << endl;
    Vector3d acc_mean = acc_sum_nav / nav_samples;
    double acc_norm = acc_mean.norm();
    
    cout << "  Local Gravity (g): " << local_g << " m/s2" << endl;
    cout << "  Meas Mean Abs(f):  " << acc_norm << " m/s2" << endl;
    cout << "  Scale Error:       " << (acc_norm - local_g) / local_g * 1e6 << " ppm" << endl;
    cout << "  Mean Acc Vector:   [" << acc_mean.transpose() << "] m/s2" << endl;

    // [New] 反推姿态 (Inferred Attitude from Accelerometer)
    // 原理：假设静止，f_b = C_n^b * [0, 0, g]^T
    // f_x = -g * sin(pitch)
    // f_y =  g * cos(pitch) * sin(roll)
    // f_z =  g * cos(pitch) * cos(roll)
    // 注意：这里的 acc_mean 包含了零偏和刻度误差，所以算出来的角度是不准的，
    // 而这个"不准"正好反映了系统当前的误差状态。
    
    double pitch_acc = asin(-acc_mean(1) / acc_norm); // 使用实测模长归一化
    double roll_acc  = asin(acc_mean(0) / acc_norm);
    
    Vector3d att_ins_end = ins.ins.att; // 惯导最终解算的姿态
    
    cout << "\n[Diagnostic] Attitude Comparison (deg):" << endl;
    cout << "  Type        Pitch       Roll" << endl;
    cout << "  --------------------------------" << endl;
    cout << "  INS (End):  " << fixed << setprecision(5) << att_ins_end(0)*glv.deg << "     " << att_ins_end(1)*glv.deg << endl;
    cout << "  Acc (Mean): " << pitch_acc*glv.deg << "     " << roll_acc*glv.deg << endl;
    cout << "  Diff:       " << (att_ins_end(0) - pitch_acc)*glv.deg << "     " << (att_ins_end(1) - roll_acc)*glv.deg << endl;

    return 0;
}