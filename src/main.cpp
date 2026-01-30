/**
 * @file main.cpp
 * @brief 最终决战：6小时稳定对准 + 24小时原子反馈导航
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
    // 1. 基础配置
    double ts = 1.0 / 400.0;
    GLV glv;
    Vector3d pos_ref(32.0286 * glv.deg, 118.8533 * glv.deg, 17.0);
    Earth eth(glv); 
    eth.update(pos_ref, Vector3d::Zero());
    double local_g = eth.gn.norm(); 

    // 文件列表 (5个分片)
    std::vector<std::string> file_list = {
        "../fog_part1.csv", 
        "../fog_part2.csv",
        "../fog_part3.csv",
        "../fog_part4.csv",
        "../fog_part5.csv"
    };

    // =========================================================
    // Phase 1: 完美复刻刚才的成功 (加载 fog_part1)
    // =========================================================
    cout << "[Phase 1] Loading Part 1: " << file_list[0] << endl;
    auto buffer_part1 = LoadIMUData(file_list[0], ts, local_g);
    
    // 初始化引擎
    SinsEngine sinsegine(ts);
    sinsegine.Sins_Init(pos_ref, Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero() ,Vector3d::Zero());

    HybridAlignConfig cfg;
    cfg.t_coarse = 60;
    // 直接用 buffer 大小，跑满 6小时，复刻刚才的成功
    cfg.t_fine = buffer_part1.size() * ts; 
    cfg.eb_sigma_allan = 0.003;
    cfg.db_sigma_allan = 50.0;
    cfg.verbose = true;
    
    cout << "Running Alignment on full Part 1 (" << cfg.t_fine << " s)..." << endl;
    
    // 这一步刚才跑通了，绝对稳
    auto align_res = sinsegine.Run_HybridAlign(buffer_part1, cfg);

    if (!align_res.valid) { 
        cerr << "Alignment failed! (Unexpected)" << endl; 
        return -1; 
    }
    
    cout << "=== ALIGNMENT SUCCESS ===" << endl;
    cout << fixed << setprecision(8);
    cout << "Att (rad):  " << align_res.att.transpose() << endl;
    cout << "Att (deg):  " << align_res.att.transpose() / glv.deg << endl;
    cout << "Gyro Bias:  " << align_res.eb.transpose() / glv.deg * 3600.0 << " deg/h" << endl;
    cout << "Acc Bias :  " << align_res.db.transpose() / 9.8 * 1e6 << " ug" << endl;

    // 释放 Part 1 内存 (因为它只用于对准，任务已完成)
    vector<IMUData>().swap(buffer_part1); 

    // =========================================================
    // Phase 2: 初始化系统 (原子 & INS)
    // =========================================================
    AtomicGyroSimulator atom(pos_ref, glv);
    atom.Init(align_res.att); // 锁定当前算出的姿态

    sinsegine.eth.update(pos_ref, Vector3d::Zero());
    // 使用对准结果重置 INS
    sinsegine.ins = INSState(align_res.att, Vector3d::Zero(), pos_ref, ts, sinsegine.eth);
    sinsegine.ins.set_bias(align_res.eb, align_res.db); 
    
    // =========================================================
    // Phase 3: 接力导航 (加载 Part 2 ~ Part 5)
    // =========================================================
    cout << "\n[Phase 3] Starting Navigation (Remaining 24h)..." << endl;
    
    ofstream log("nav_30h_final.csv");
    log << "time,lat_err_m,lon_err_m,h_err_m,drift_m,vn,ve,vd,roll,pitch,yaw\n";

    double max_drift = 0;
    size_t total_steps = 0;
    
    // 原子反馈参数
    double T_cycle = 2.0;
    int samples_per_cycle = static_cast<int>(T_cycle / ts);
    Vector3d fog_sum = Vector3d::Zero();
    int sample_count = 0;

    auto ProcessEpoch = [&](const IMUData& epoch) {
        // 1. 累积原子反馈
        fog_sum += epoch.wm / ts; 
        sample_count++;

        // 2. 惯导解算
        sinsegine.Step_Nav(epoch.wm, epoch.vm);

        // 3. 反馈执行
        if (sample_count >= samples_per_cycle) {
            Vector3d fog_mean = fog_sum / sample_count;
            Vector3d atom_meas = atom.Measure(); 
            Vector3d bias_new = fog_mean - atom_meas;
            // 核心：保留对准算出的加计零偏，只修陀螺
            sinsegine.ins.set_bias(bias_new, align_res.db); 
            fog_sum.setZero();
            sample_count = 0;
        }

        // 4. 日志记录 (1Hz)
        if (total_steps % 400 == 0) { 
            double nav_time = total_steps * ts; 
            double lat_err = (sinsegine.ins.pos(0) - pos_ref(0)) * glv.Re;
            double lon_err = (sinsegine.ins.pos(1) - pos_ref(1)) * glv.Re * cos(pos_ref(0));
            double h_err = sinsegine.ins.pos(2) - pos_ref(2);
            double drift = sqrt(lat_err*lat_err + lon_err*lon_err);
            if(drift > max_drift) max_drift = drift;

            Vector3d att_deg = sinsegine.ins.att / glv.deg; 
            Vector3d vn = sinsegine.ins.vn;
            
            log << fixed << setprecision(6) 
                << nav_time << "," 
                << lat_err << "," << lon_err << "," << h_err << "," << drift << ","
                << vn(0) << "," << vn(1) << "," << vn(2) << ","     
                << att_deg(0) << "," << att_deg(1) << "," << att_deg(2) << "\n";
            
            // 每小时打印一次进度
            if (total_steps > 0 && total_steps % (3600 * 400) == 0) {
                 cout << "Nav Time: " << nav_time / 3600.0 << " h | Drift: " << max_drift << " m | Yaw: " << att_deg(2) << endl;
            }
        }
        total_steps++;
    };

    // 使用 Loader 加载 Part 2 ~ Part 5 (从第2个文件开始)
    cout << "Switching to ChainedLoader for Part 2-5..." << endl;
    std::vector<std::string> rest_files(file_list.begin() + 1, file_list.end());
    IMUChainedLoader nav_loader(rest_files, ts, local_g);
    IMUData epoch;

    while (nav_loader.Next(epoch)) {
        ProcessEpoch(epoch);
    }

    log.close();
    cout << "--------------------------------" << endl;
    cout << "Total Nav Time: " << total_steps * ts / 3600.0 << " h" << endl;
    cout << "Max Drift:      " << max_drift << " m" << endl;
    cout << "Final Drift:    " << max_drift / 1852.0 << " nm" << endl;

    return 0;
}