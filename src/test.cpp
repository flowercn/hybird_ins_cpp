/**
 * @file test.cpp
 * @brief 纯惯导测试 - 使用封装后的简洁接口
 * 
 * 接口说明：
 *   - LoadIMUData(file, ts, g): 加载 IMU 数据
 *   - SinsEngine::Run_HybridAlign(data, cfg): 一键混合对准
 *   - SinsEngine::Step_Nav(wm, vm): 单步导航
 * 
 * 混合对准原理：
 *   1. 粗对准：解析法求初始姿态
 *   2. KF 精对准：估计精确姿态（忽略其零偏估计）
 *   3. 几何均值法：直接从陀螺数据估计零偏
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include "sins_engine.h"
#include "support.h"

using namespace std;
using namespace Eigen;

int main() {
    // // ========== 参数配置 ==========
    // double ts = 1.0 / 400.0;
    // GLV glv;
    
    // // 初始位置
    // Vector3d pos(32.0286 * glv.deg, 118.8533 * glv.deg, 17.0);
    
    // // 当地重力
    // Earth eth(glv);
    // eth.update(pos, Vector3d::Zero());
    // double local_g = eth.gn.norm();
    
    // // ========== 加载数据 ==========
    // cout << "[1] Loading data..." << endl;
    // auto imu_data = LoadIMUData("../fog3h.csv", ts, local_g);
    // cout << "    " << imu_data.size() * ts / 60 << " min loaded" << endl;
    
    // // ========== 初始化 + 对准 ==========
    // cout << "\n[2] Initializing and aligning..." << endl;
    
    // SinsEngine ins(ts);
    // ins.Sins_Init(pos, Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero());
    
    // HybridAlignConfig cfg;
    // cfg.t_coarse = 60;          // 60s 粗对准
    // cfg.t_fine = 3600;          // 1h 精对准
    // cfg.eb_sigma_allan = 0.003; // Allan 零偏稳定性 (deg/h)
    // cfg.verbose = true;
    
    // auto result = ins.Run_HybridAlign(imu_data, cfg);
    
    // if (!result.valid) {
    //     cerr << "Alignment failed!" << endl;
    //     return -1;
    // }
    
    // // ========== 导航 ==========
    // cout << "\n[3] Navigation..." << endl;
    
    // // 重置状态
    // ins.eth.update(pos, Vector3d::Zero());
    // ins.ins = INSState(result.att, Vector3d::Zero(), pos, ts, ins.eth);
    // ins.ins.set_bias(result.eb, result.db);
    
    // // 打开日志文件
    // ofstream log_file("nav_log.csv");
    // log_file << "time_s,lat_deg,lon_deg,alt_m,pitch_deg,roll_deg,yaw_deg,drift_m\n";
    
    // // 导航
    // size_t nav_start = static_cast<size_t>(cfg.t_fine / ts);
    // double max_drift = 0;
    
    // cout << "    Time(min)  Drift(m)" << endl;
    // int step = (imu_data.size() - nav_start) / 10;
    // int log_step = static_cast<int>(1.0 / ts);  // 每秒记录一次
    
    // for (size_t i = nav_start; i < imu_data.size(); ++i) {
    //     ins.Step_Nav(imu_data[i].wm, imu_data[i].vm);
        
    //     double dN = (ins.ins.pos(0) - pos(0)) * glv.Re;
    //     double dE = (ins.ins.pos(1) - pos(1)) * cos(pos(0)) * glv.Re;
    //     double drift = sqrt(dN*dN + dE*dE);
    //     if (drift > max_drift) max_drift = drift;
        
    //     // 保存日志 (每秒)
    //     if ((i - nav_start) % log_step == 0) {
    //         double t = (i - nav_start) * ts;
    //         Vector3d att = INSMath::m2att(ins.ins.Cnb);
    //         log_file << fixed << setprecision(6)
    //                  << t << ","
    //                  << ins.ins.pos(0) / glv.deg << ","
    //                  << ins.ins.pos(1) / glv.deg << ","
    //                  << ins.ins.pos(2) << ","
    //                  << att(0) / glv.deg << ","
    //                  << att(1) / glv.deg << ","
    //                  << att(2) / glv.deg << ","
    //                  << drift << "\n";
    //     }
        
    //     if (step > 0 && (i - nav_start) % step == 0 && i > nav_start) {
    //         cout << fixed << setw(12) << setprecision(1) << ((i - nav_start) * ts / 60)
    //              << setw(10) << drift << endl;
    //     }
    // }
    
    // log_file.close();
    // cout << "    Log saved to: nav_log.csv" << endl;
    
    // // ========== 结果 ==========
    // double nav_time = (imu_data.size() - nav_start) * ts;
    // double drift_rate = max_drift / (nav_time / 3600) / 1852;
    
    // cout << "\n[4] Result:" << endl;
    // cout << "    Nav time:    " << nav_time / 60 << " min" << endl;
    // cout << "    Max drift:   " << fixed << setprecision(1) << max_drift << " m" << endl;
    // cout << "    Drift rate:  " << setprecision(3) << drift_rate << " nm/h" << endl;
    
    // return 0;
}
