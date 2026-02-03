#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "sins_engine.h"
#include "cai_sim.h"
#include "support.h"

using namespace std;
using namespace Eigen;

struct ExpConfig {
    string name;        // 实验名称
    double vrw_ug;      // 原子加计白噪声
    double bias_ug;     // 原子加计常值零偏
    double gain;        // 反馈增益
    string filename;    // 输出文件名
};


int main() {
    // 1. 基础配置
    double ts = 1.0 / 400.0;
    GLV glv;
    Vector3d pos_ref(32.0286 * glv.deg, 118.8533 * glv.deg, 17.0);
    Earth eth(glv); 
    eth.update(pos_ref, Vector3d::Zero());
    double local_g = eth.gn.norm(); 

    std::vector<std::string> file_list = {
        "../fog_part1.csv", "../fog_part2.csv", "../fog_part3.csv", 
        "../fog_part4.csv", "../fog_part5.csv"
    };

    // =========================================================
    // Phase 0: 纯原子导航基准测试
    // =========================================================
    {
        Vector3d att_true = Vector3d(0.0436479, 0.332121, 0.553816) * glv.deg;
        AtomicGyroSimulator atom_pure_g(pos_ref, glv);
        atom_pure_g.Init(att_true);

        CAIParams acc_pure_p;
        acc_pure_p.bias_ug = 0.0;
        acc_pure_p.vrw_ug  = 0.05; 
        AtomicAccSimulator atom_pure_a(pos_ref, glv, acc_pure_p);
        atom_pure_a.Init(att_true);

        // 运行测试 (注意：这里去掉了 exit(0)，让它跑完继续)
        Run_Pure_Atomic_Nav(pos_ref, att_true, atom_pure_g, atom_pure_a, ts);
    }

    // =========================================================
    // Phase 1: 10分钟对准 (使用 5h50m ~ 6h00m 数据)
    // =========================================================
    cout << "[Phase 1] Loading Part 1: " << file_list[0] << endl;
    auto buffer_part1 = LoadIMUData(file_list[0], ts, local_g);
    
    // --- 提取最后 10 分钟数据 ---
    double start_time_sec = 5.0 * 3600.0 + 50.0 * 60.0; // 5h50m = 21000s
    size_t start_idx = static_cast<size_t>(start_time_sec / ts);
    
    cout << ">>> Slicing data from " << start_time_sec << "s to end (Last 10 mins)..." << endl;
    std::vector<IMUData> buffer_align(buffer_part1.begin() + start_idx, buffer_part1.end());
    cout << ">>> Alignment Data Size: " << buffer_align.size() << " epochs (" << buffer_align.size() * ts << " s)" << endl;


    // 手动初始化 (绕过 Sins_Init 的潜在 Bug)
    SinsEngine sinsegine(ts);
    sinsegine.eth = Earth(glv);
    sinsegine.eth.update(pos_ref, Vector3d::Zero());
    sinsegine.ins = INSState(Vector3d::Zero(), Vector3d::Zero(), pos_ref, ts, sinsegine.eth);
    sinsegine.res_init.pos = pos_ref; // 确保 res_init 也有值

    HybridAlignConfig cfg;
    cfg.t_coarse = 60;
    cfg.t_fine = buffer_align.size() * ts; 
    cfg.eb_sigma_allan = 0.003;
    cfg.db_sigma_allan = 50.0;
    cfg.verbose = true;
    
    cout << "Running Alignment on full Part 1 (" << cfg.t_fine << " s)..." << endl;
    auto align_res = sinsegine.Run_HybridAlign(buffer_align, cfg);
    if (!align_res.valid) { cerr << "Alignment failed!" << endl; return -1; }
    
    cout << "=== ALIGNMENT SUCCESS ===" << endl;
    cout << "Att (deg):  " << align_res.att.transpose() / glv.deg << endl;
    cout << "Acc Bias :  " << align_res.db.transpose() / glv.ug << " ug" << endl;

    // 释放大内存
    vector<IMUData>().swap(buffer_part1); 
    vector<IMUData>().swap(buffer_align);
    
    // =========================================================
    // Phase 2: 初始化原子系统 (参数配置：10ug, Gain 1.0)
    // =========================================================

    // 原子陀螺
    AtomicGyroSimulator atom(pos_ref, glv);
    atom.Init(align_res.att); 

    // 原子加计 
    CAIParams acc_params;
    acc_params.bias_ug = 0.0; 
    acc_params.vrw_ug  = 0.05; // 关键：低噪声
    AtomicAccSimulator atom_acc(pos_ref, glv, acc_params);
    atom_acc.Init(align_res.att); 

    vector<IMUData>().swap(buffer_part1); // 释放内存

    // =========================================================
    // Phase 3: 24小时导航
    // =========================================================
    double gain = 1.0; // 【请确保这里的 gain 和下面用的一致】
    double T_cycle = 2.0;

    cout << "\n==================================================" << endl;
    cout << "           NAVIGATION CONFIGURATION               " << endl;
    cout << "==================================================" << endl;
    cout << " [Sensor] Atomic Acc Bias : " << acc_params.bias_ug << " ug" << endl;
    cout << " [Sensor] Atomic Acc VRW  : " << acc_params.vrw_ug  << " ug" << endl;
    cout << " [Algo  ] Feedback Gain   : " << gain << endl;
    cout << " [Algo  ] Atomic Cycle    : " << T_cycle << " s" << endl;
    cout << " [Align ] Alignment Time  : " << cfg.t_fine << " s" << endl;
    cout << " [Output] CSV File        : nav_30h_final.csv" << endl;
    cout << "==================================================\n" << endl;

    cout << "[Phase 3] Starting Navigation (Remaining 24h)..." << endl;

    // 1. 手动初始化 INS
    sinsegine.eth.update(pos_ref, Vector3d::Zero());
    sinsegine.ins = INSState(align_res.att, Vector3d::Zero(), pos_ref, ts, sinsegine.eth);
    sinsegine.ins.set_bias(align_res.eb, align_res.db);

    ofstream log("nav_30h_final.csv");
    // 请勿修改表头
    log << "time,lat_err,lon_err,h_err,drift,vn,ve,vu,roll,pitch,yaw,"
        << "eb_x,eb_y,eb_z,db_x,db_y,db_z\n";

    double max_drift = 0;
    size_t total_steps = 0;
    
    int samples_per_cycle = static_cast<int>(T_cycle / ts);
    std::vector<IMUData> chunk_buffer;
    chunk_buffer.reserve(samples_per_cycle + 10);
    
    Vector3d fog_gyro_sum = Vector3d::Zero();
    Vector3d fog_vel_sum  = Vector3d::Zero();
    int sample_count = 0;

    

    auto ProcessEpoch = [&](const IMUData& epoch) {
        chunk_buffer.push_back(epoch);
        fog_gyro_sum += epoch.wm / ts; 
        fog_vel_sum  += epoch.vm;     
        sample_count++;

        if (sample_count >= samples_per_cycle) {
            Vector3d fog_gyro_mean = fog_gyro_sum / sample_count;
            Vector3d fog_acc_mean  = fog_vel_sum / (sample_count * ts);
            Vector3d atom_gyro_meas = atom.Measure(); 
            Vector3d atom_acc_meas  = atom_acc.Measure(); 

            Vector3d eb_obs = fog_gyro_mean - atom_gyro_meas;
            Vector3d db_obs = fog_acc_mean - atom_acc_meas;

            Vector3d eb_new = sinsegine.ins.eb * (1.0 - gain) + eb_obs * gain;
            Vector3d db_new = sinsegine.ins.db * (1.0 - gain) + db_obs * gain;
            sinsegine.ins.set_bias(eb_new, db_new);

            for (const auto& buffered_epoch : chunk_buffer) {
                sinsegine.Step_Nav(buffered_epoch.wm, buffered_epoch.vm);
                
                // 【阻尼】锁定天向 (ENU系下 2=Up)
                sinsegine.ins.vn(2) = 0.0; 

                if (total_steps % 400 == 0) { 
                    double nav_time = total_steps * ts; 
                    double lat_err = (sinsegine.ins.pos(0) - pos_ref(0)) * glv.Re;
                    double lon_err = (sinsegine.ins.pos(1) - pos_ref(1)) * glv.Re * cos(pos_ref(0));
                    double h_err = sinsegine.ins.pos(2) - pos_ref(2);
                    double drift = sqrt(lat_err*lat_err + lon_err*lon_err);
                    if(drift > max_drift) max_drift = drift;
                    
                    // 请勿删除
                    log << fixed << setprecision(6) 
                        << nav_time << "," 
                        << lat_err << "," 
                        << lon_err << "," 
                        << h_err << "," 
                        << drift << ","
                        << sinsegine.ins.vn(0) << ","  // Vn (East)
                        << sinsegine.ins.vn(1) << ","  // Ve (North)
                        << sinsegine.ins.vn(2) << ","  // Vu (Up)
                        << sinsegine.ins.att(0)/glv.deg << "," // Pitch
                        << sinsegine.ins.att(1)/glv.deg << "," // Roll
                        << sinsegine.ins.att(2)/glv.deg << "," // Yaw
                        
                        // ---陀螺零偏 (deg/h) ---
                        << sinsegine.ins.eb(0)/glv.deg * 3600.0 << "," // eb_x
                        << sinsegine.ins.eb(1)/glv.deg * 3600.0 << "," // eb_y
                        << sinsegine.ins.eb(2)/glv.deg * 3600.0 << "," // eb_z
                        
                        // --- 加计零偏 (ug) ---
                        << sinsegine.ins.db(0)/glv.ug << ","       // db_x
                        << sinsegine.ins.db(1)/glv.ug << ","       // db_y
                        << sinsegine.ins.db(2)/glv.ug << "\n";     // db_z
                    
                    if (total_steps > 0 && total_steps % (3600 * 400) == 0) {
                         cout << "Time: " << nav_time / 3600.0 << " h | Drift: " << drift 
                              << " m | Max: " << max_drift << " m" << endl;
                    }
                }
                total_steps++;
            }
            fog_gyro_sum.setZero(); fog_vel_sum.setZero(); sample_count = 0; chunk_buffer.clear();
        }
    };

    cout << "Switching to ChainedLoader..." << endl;
    std::vector<std::string> rest_files(file_list.begin() + 1, file_list.end());
    IMUChainedLoader nav_loader(rest_files, ts, local_g);
    IMUData epoch;
    while (nav_loader.Next(epoch)) ProcessEpoch(epoch);

    log.close();
    cout << "Final Drift: " << max_drift << " m" << endl;
    return 0;
}