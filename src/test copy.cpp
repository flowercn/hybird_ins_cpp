#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include "sins_engine.h"
#include "cai_sim.h"
#include "support.h" // 包含 ExpConfig 等基础定义

using namespace std;
using namespace Eigen;
using namespace cai;

// 定义死区时间实验配置
struct DeadTimeConfig {
    string name;        // 实验组名称
    double t_active;    // 原子陀螺有效工作时间
    double t_dead;      // 盲区时间
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
    // Phase 1: 10分钟对准 (共用)
    // =========================================================
    cout << "\n[Phase 1] Loading Alignment Data: " << file_list[0] << endl;
    auto buffer_part1 = LoadIMUData(file_list[0], ts, local_g);
    
    // 提取最后 10 分钟
    double start_time_sec = 5.0 * 3600.0 + 50.0 * 60.0; 
    size_t start_idx = static_cast<size_t>(start_time_sec / ts);
    std::vector<IMUData> buffer_align(buffer_part1.begin() + start_idx, buffer_part1.end());
    
    SinsEngine align_engine(ts);
    align_engine.eth = Earth(glv);
    align_engine.eth.update(pos_ref, Vector3d::Zero());
    align_engine.ins = INSState(Vector3d::Zero(), Vector3d::Zero(), pos_ref, ts, align_engine.eth);
    align_engine.res_init.pos = pos_ref;

    HybridAlignConfig cfg;
    cfg.t_coarse = 60;
    cfg.t_fine = buffer_align.size() * ts; 
    cfg.eb_sigma_allan = 0.003;
    cfg.db_sigma_allan = 50.0;
    cfg.verbose = false; 
    
    cout << "Running Alignment (" << cfg.t_fine << " s)..." << endl;
    auto align_res = align_engine.Run_HybridAlign(buffer_align, cfg);
    if (!align_res.valid) { cerr << "Alignment failed!" << endl; return -1; }
    
    cout << "=== ALIGNMENT SUCCESS ===" << endl;
    cout << "Att (deg):  " << align_res.att.transpose() / glv.deg << endl;
    cout << "Acc Bias :  " << align_res.db.transpose() / glv.ug << " ug" << endl;

    vector<IMUData>().swap(buffer_part1); 
    vector<IMUData>().swap(buffer_align);
    
    // =========================================================
    // Phase 3: 6组死区时间对比实验 (手动回溯版 + 盲区残差记录)
    // =========================================================

    vector<DeadTimeConfig> experiments = {
        {"Group_2.0s_Full", 2.0,      0.0,      "nav_dead_2.0s.csv"}, 
        //{"Group_1.8s",      1.8,      0.2,      "nav_dead_1.8s.csv"},
        {"Group_1.6s",      1.6,      0.4,      "nav_dead_1.6s.csv"}, 
       // {"Group_1.4s",      1.4,      0.6,      "nav_dead_1.4s.csv"},
        {"Group_1.2s",      1.2,      0.8,      "nav_dead_1.2s.csv"},
        //{"Group_1.0s",      1.0,      1.0,      "nav_dead_1.0s.csv"}  
    };

    std::vector<std::string> rest_files(file_list.begin() + 1, file_list.end());
    cout << "\n>>> Starting Dead-Time Analysis (" << experiments.size() << " Scenarios)..." << endl;

    for (const auto& exp : experiments) {
        cout << "\n--------------------------------------------------" << endl;
        cout << " Running Scenario: " << exp.name << endl;
        cout << " Config: Active=" << exp.t_active << "s, Dead=" << exp.t_dead << "s -> " << exp.filename << endl;
        cout << "--------------------------------------------------" << endl;

        // 1. 初始化传感器
        AtomicGyroSimulator atom(pos_ref, glv);
        atom.Init(align_res.att); 
        
        // 2. 初始化惯导引擎
        SinsEngine sinsegine(ts);
        sinsegine.eth = Earth(glv);
        sinsegine.eth.update(pos_ref, Vector3d::Zero());
        sinsegine.ins = INSState(align_res.att, Vector3d::Zero(), pos_ref, ts, sinsegine.eth);
        sinsegine.ins.set_bias(align_res.eb, align_res.db);

        // 3. 日志配置
        IMUChainedLoader nav_loader(rest_files, ts, local_g);
        ofstream log(exp.filename);
        // [修改] 增加残差列 res_int_x, res_int_y, res_int_z
        log << "time,lat_err,lon_err,h_err,drift,vn,ve,vu,roll,pitch,yaw,"
            << "eb_x,eb_y,eb_z,db_x,db_y,db_z,"
            << "res_int_x,res_int_y,res_int_z\n"; 

        // 4. 周期参数
        double T_cycle = 2.0; 
        int total_samples_per_cycle = static_cast<int>(T_cycle / ts); 
        int active_samples = static_cast<int>(exp.t_active / ts);     
        
        std::vector<IMUData> chunk_buffer;
        chunk_buffer.reserve(total_samples_per_cycle + 10);
        
        Vector3d fog_gyro_active_sum = Vector3d::Zero(); 
        
        // 【核心变量】盲区残差累加器 (Sum of Raw - Earth)
        Vector3d dead_residual_sum = Vector3d::Zero(); 
        
        int sample_count = 0;
        double max_drift = 0;
        size_t total_steps = 0;
        double gain = 1.0; 

        auto ProcessEpoch = [&](const IMUData& epoch) {
            // A. 数据缓存
            chunk_buffer.push_back(epoch);
            sample_count++;

            // B. 有条件累加
            if (sample_count <= active_samples) {
                // 有效期：累加 FOG 用于估算 Bias
                fog_gyro_active_sum += epoch.wm / ts; 
            } else {
                // 【核心修改】盲区：计算残差 (光纤原始输出 - 地球自转分量)
                
                // 1. 计算当前时刻(近似)的地球自转在载体系的分量
                //    使用周期起始时的姿态 Cnb 进行投影
                Vector3d earth_inc_body = sinsegine.ins.Cnb.transpose() * sinsegine.eth.wien * ts;
                
                // 2. 累加残差 (Res = Wm - Earth)
                //    这里剩下的就是：Bias * dt + Noise
                dead_residual_sum += (epoch.wm - earth_inc_body);
            }

            // D. 周期结束 -> 计算与回溯
            if (sample_count >= total_samples_per_cycle) {
                
                // 1. 估算 Active 区 Bias
                Vector3d fog_gyro_mean = Vector3d::Zero();
                if (active_samples > 0) fog_gyro_mean = fog_gyro_active_sum / active_samples;
                Vector3d atom_gyro_meas = atom.Measure(); 
                Vector3d eb_obs = fog_gyro_mean - atom_gyro_meas;
                Vector3d eb_new = sinsegine.ins.eb * (1.0 - gain) + eb_obs * gain;
                sinsegine.ins.set_bias(eb_new, sinsegine.ins.db); 

                // 2. 锁定本周期的盲区残差，用于后续写入日志
                Vector3d current_res_int = dead_residual_sum;

                // 3. [回溯] 全周期重算
                for (size_t i = 0; i < chunk_buffer.size(); ++i) {
                    const auto& buffered_epoch = chunk_buffer[i];
                    sinsegine.Step_Nav(buffered_epoch.wm, buffered_epoch.vm);
                    sinsegine.ins.vn(2) = 0.0; // 阻尼

                    if (total_steps % 400 == 0) { 
                        double nav_time = total_steps * ts; 
                        double lat_err = (sinsegine.ins.pos(0) - pos_ref(0)) * glv.Re;
                        double lon_err = (sinsegine.ins.pos(1) - pos_ref(1)) * glv.Re * cos(pos_ref(0));
                        double h_err = sinsegine.ins.pos(2) - pos_ref(2);
                        double drift = sqrt(lat_err*lat_err + lon_err*lon_err);
                        if(drift > max_drift) max_drift = drift;

                        log << fixed << setprecision(15) 
                            << nav_time << "," << lat_err << "," << lon_err << "," << h_err << "," << drift << ","
                            << sinsegine.ins.vn(0) << "," << sinsegine.ins.vn(1) << "," << sinsegine.ins.vn(2) << "," 
                            << sinsegine.ins.att(0)/glv.deg << "," << sinsegine.ins.att(1)/glv.deg << "," << sinsegine.ins.att(2)/glv.deg << ","
                            << sinsegine.ins.eb(0)/glv.deg*3600 << "," << sinsegine.ins.eb(1)/glv.deg*3600 << "," << sinsegine.ins.eb(2)/glv.deg*3600 << ","
                            << sinsegine.ins.db(0)/glv.ug << "," << sinsegine.ins.db(1)/glv.ug << "," << sinsegine.ins.db(2)/glv.ug << ","
                            // 输出盲区残差 (rad)
                            << current_res_int(0) << "," << current_res_int(1) << "," << current_res_int(2) << "\n";
                        
                        if (total_steps > 0 && total_steps % (3600 * 400) == 0 && i == chunk_buffer.size()-1) {
                             cout << "   Time: " << nav_time / 3600.0 << " h | Drift: " << drift 
                                  << " m | Max: " << max_drift << " m" << endl;
                        }
                    }
                    total_steps++;
                }

                // 清空状态
                fog_gyro_active_sum.setZero();
                dead_residual_sum.setZero(); // 重置残差累加器
                sample_count = 0; 
                chunk_buffer.clear();
            }
        };

        IMUData epoch;
        while (nav_loader.Next(epoch)) ProcessEpoch(epoch);
        log.close();
        cout << " [Done] Final Drift: " << max_drift << " m" << endl;
    } 

    cout << "\nAll Dead-Time scenarios completed successfully!" << endl;
    return 0;
}