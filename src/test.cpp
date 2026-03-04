#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include "sins_engine.h"
#include "cai_sim.h"
#include "support.h" 

using namespace std;
using namespace Eigen;
using namespace cai;

// ==========================================
// 单步运行函数
// ==========================================
void RunSingleExperiment(const ExperimentConfig& exp, const SimulationContext& ctx) {
    cout << "\n--------------------------------------------------" << endl;
    cout << " Running Scenario: " << exp.name << endl;
    cout << " Config: Active=" << exp.t_active << "s, Dead=" << exp.t_dead << "s -> " << exp.output_file << endl;
    cout << "--------------------------------------------------" << endl;

    // A. 初始化传感器 (原子陀螺 & 原子加计)
    AtomicGyroSimulator atom(ctx.pos_ref, ctx.glv);
    atom.Init(ctx.att_align); 
    
    CAIParams acc_params;
    acc_params.bias_ug = 0.0; 
    acc_params.vrw_ug  = 0.05; 
    AtomicAccSimulator atom_acc(ctx.pos_ref, ctx.glv, acc_params);
    atom_acc.Init(ctx.att_align); 
    
    // B. 初始化惯导引擎
    SinsEngine sinsegine(ctx.ts);
    sinsegine.eth = ctx.eth; // 直接拷贝预设好的 Earth 对象
    sinsegine.ins = INSState(ctx.att_align, Vector3d::Zero(), ctx.pos_ref, ctx.ts, sinsegine.eth);
    sinsegine.ins.set_bias(ctx.eb_align, ctx.db_align);

    // C. 准备数据加载器 (跳过第一个对准文件)
    std::vector<std::string> rest_files(ctx.file_list.begin() + 1, ctx.file_list.end());
    IMUChainedLoader nav_loader(rest_files, ctx.ts, ctx.local_g);

    // D. 准备日志
    ofstream log(exp.output_file);
    log << "time,lat_err,lon_err,h_err,drift,vn,ve,vu,roll,pitch,yaw,"
        << "eb_x,eb_y,eb_z,db_x,db_y,db_z,"
        << "res_gyro_x,res_gyro_y,res_gyro_z,"
        << "res_acc_x,res_acc_y,res_acc_z\n"; 

    // E. 循环变量
    double T_cycle = 2.0; 
    int total_samples_per_cycle = static_cast<int>(T_cycle / ctx.ts); 
    int active_samples = static_cast<int>(exp.t_active / ctx.ts);     
    
    std::vector<IMUData> chunk_buffer;
    chunk_buffer.reserve(total_samples_per_cycle + 10);
    
    Vector3d fog_gyro_active_sum = Vector3d::Zero(); 
    Vector3d fog_acc_active_sum  = Vector3d::Zero(); 
    Vector3d gyro_dead_raw_sum = Vector3d::Zero(); 
    Vector3d acc_dead_raw_sum  = Vector3d::Zero(); 
    
    int sample_count = 0;
    int dead_count = 0;
    double max_drift = 0;
    size_t total_steps = 0;

    // F. 处理每一帧数据的 Lambda
    auto ProcessEpoch = [&](const IMUData& epoch) {
        chunk_buffer.push_back(epoch);
        sample_count++;

        if (sample_count <= active_samples) {
            // [Active]
            fog_gyro_active_sum += epoch.wm / ctx.ts; 
            fog_acc_active_sum  += epoch.vm / ctx.ts; 
        } else {
            // [Dead] 关键修改点：使用 ctx.Cnb_align 作为固定基准
            Vector3d earth_inc_body = ctx.Cnb_align.transpose() * ctx.eth.wien * ctx.ts;
            gyro_dead_raw_sum += (epoch.wm - earth_inc_body);

            Vector3d force_inc_body = ctx.Cnb_align.transpose() * (-ctx.eth.gn) * ctx.ts;
            acc_dead_raw_sum += (epoch.vm - force_inc_body);

            dead_count++;
        }

        // 周期结算
        if (sample_count >= total_samples_per_cycle) {
            // 1. 估算 Bias
            Vector3d fog_gyro_mean = Vector3d::Zero();
            Vector3d fog_acc_mean  = Vector3d::Zero();
            if (active_samples > 0) {
                fog_gyro_mean = fog_gyro_active_sum / active_samples;
                fog_acc_mean  = fog_acc_active_sum  / active_samples;
            }

            Vector3d atom_gyro_meas = atom.Measure(); 
            Vector3d eb_new = fog_gyro_mean - atom_gyro_meas;

            Vector3d db_new = sinsegine.ins.db;
            if (exp.use_atomic_acc && active_samples > 0) {
                Vector3d atom_acc_meas  = atom_acc.Measure(); 
                db_new = fog_acc_mean  - atom_acc_meas;
            }

            sinsegine.ins.set_bias(eb_new, db_new);

            // 2. 计算残差
            double t_dead_real = dead_count * ctx.ts;
            Vector3d final_gyro_res = gyro_dead_raw_sum - eb_new * t_dead_real;
            Vector3d final_acc_res  = acc_dead_raw_sum  - db_new * t_dead_real;

            // 3. 回溯
            for (size_t i = 0; i < chunk_buffer.size(); ++i) {
                const auto& buf_epoch = chunk_buffer[i];
                sinsegine.Step_Nav(buf_epoch.wm, buf_epoch.vm);
                sinsegine.ins.vn(2) = 0.0; // 阻尼

                if (total_steps % 400 == 0) { 
                    double nav_time = total_steps * ctx.ts; 
                    
                    // --- 完整的误差计算逻辑 (无删减) ---
                    double lat_err = (sinsegine.ins.pos(0) - ctx.pos_ref(0)) * ctx.glv.Re;
                    double lon_err = (sinsegine.ins.pos(1) - ctx.pos_ref(1)) * ctx.glv.Re * cos(ctx.pos_ref(0));
                    double h_err = sinsegine.ins.pos(2) - ctx.pos_ref(2);
                    double drift = sqrt(lat_err*lat_err + lon_err*lon_err);
                    if(drift > max_drift) max_drift = drift;

                    log << fixed << setprecision(15) 
                        << nav_time << "," << lat_err << "," << lon_err << "," << h_err << "," << drift << ","
                        << sinsegine.ins.vn(0) << "," << sinsegine.ins.vn(1) << "," << sinsegine.ins.vn(2) << "," 
                        << sinsegine.ins.att(0)/ctx.glv.deg << "," << sinsegine.ins.att(1)/ctx.glv.deg << "," << sinsegine.ins.att(2)/ctx.glv.deg << ","
                        << sinsegine.ins.eb(0)/ctx.glv.deg*3600 << "," << sinsegine.ins.eb(1)/ctx.glv.deg*3600 << "," << sinsegine.ins.eb(2)/ctx.glv.deg*3600 << ","
                        << sinsegine.ins.db(0)/ctx.glv.ug << "," << sinsegine.ins.db(1)/ctx.glv.ug << "," << sinsegine.ins.db(2)/ctx.glv.ug << ","
                        << final_gyro_res(0) << "," << final_gyro_res(1) << "," << final_gyro_res(2) << ","
                        << final_acc_res(0)  << "," << final_acc_res(1)  << "," << final_acc_res(2) << "\n";
                    
                    // 仅每小时打印进度
                    if (total_steps > 0 && total_steps % (3600 * 400) == 0 && i == chunk_buffer.size()-1) {
                         cout << "   Time: " << nav_time / 3600.0 << " h | Drift: " << drift 
                              << " m | Max: " << max_drift << " m" << endl;
                    }
                }
                total_steps++;
            }

            // 重置
            fog_gyro_active_sum.setZero(); fog_acc_active_sum.setZero();
            gyro_dead_raw_sum.setZero();   acc_dead_raw_sum.setZero();
            dead_count = 0;
            sample_count = 0; 
            chunk_buffer.clear();
        }
    };

    IMUData epoch;
    while (nav_loader.Next(epoch)) ProcessEpoch(epoch);
    log.close();
    cout << " [Done] Final Drift: " << max_drift << " m" << endl;
}

// ==========================================
// 主函数
// ==========================================
int main() {
    // 1. 准备环境 Context
    SimulationContext ctx;
    ctx.ts = 1.0 / 400.0;
    ctx.glv = GLV(); // 确保GLV初始化
    ctx.pos_ref = Vector3d(32.0286 * ctx.glv.deg, 118.8533 * ctx.glv.deg, 17.0);
    ctx.eth = Earth(ctx.glv);
    ctx.eth.update(ctx.pos_ref, Vector3d::Zero());
    ctx.local_g = ctx.eth.gn.norm();
    
    ctx.file_list = {
        "../fog_part1.csv", "../fog_part2.csv", "../fog_part3.csv", 
        "../fog_part4.csv", "../fog_part5.csv"
    };

    // 2. 执行对准 (Phase 1)
    cout << "\n[Phase 1] Loading Alignment Data: " << ctx.file_list[0] << endl;
    auto buffer_part1 = LoadIMUData(ctx.file_list[0], ctx.ts, ctx.local_g);
    
    double start_time_sec = 5.0 * 3600.0 + 50.0 * 60.0; 
    size_t start_idx = static_cast<size_t>(start_time_sec / ctx.ts);
    std::vector<IMUData> buffer_align(buffer_part1.begin() + start_idx, buffer_part1.end());

    SinsEngine align_engine(ctx.ts);
    align_engine.eth = ctx.eth;
    align_engine.ins = INSState(Vector3d::Zero(), Vector3d::Zero(), ctx.pos_ref, ctx.ts, ctx.eth);
    align_engine.res_init.pos = ctx.pos_ref;

    HybridAlignConfig cfg;
    cfg.t_coarse = 60;
    cfg.t_fine = buffer_align.size() * ctx.ts; 
    cfg.eb_sigma_allan = 0.003;
    cfg.db_sigma_allan = 50.0;
    cfg.verbose = false; 

    cout << "Running Alignment (" << cfg.t_fine << " s)..." << endl;
    auto align_res = align_engine.Run_HybridAlign(buffer_align, cfg);
    if (!align_res.valid) { cerr << "Alignment failed!" << endl; return -1; }
    
    cout << "=== ALIGNMENT SUCCESS ===" << endl;
    cout << "Att (deg):  " << align_res.att.transpose() / ctx.glv.deg << endl;
    cout << "Acc Bias :  " << align_res.db.transpose() / ctx.glv.ug << " ug" << endl;

    // 填充 Context 中的对准结果
    ctx.att_align = align_res.att;
    ctx.eb_align  = align_res.eb;
    ctx.db_align  = align_res.db;
    // 【关键】锁定对准后的姿态矩阵
    ctx.Cnb_align = INSMath::a2mat(align_res.att); 

    // 释放内存
    vector<IMUData>().swap(buffer_part1); 
    vector<IMUData>().swap(buffer_align);

    // 3. 批量运行实验 (Phase 3)
    vector<ExperimentConfig> experiments = {
        // 第一组：仅陀螺组合 (复现舒拉震荡)
        //{"Group_2.0s_GyroOnly", 2.0, 0.0, false, "nav_dead_2.0s_gyro_only.csv"}, 
        //{"Group_1.6s_GyroOnly", 1.6, 0.4, false, "nav_dead_1.6s_gyro_only.csv"}, 
        //{"Group_1.2s_GyroOnly", 1.2, 0.8, false, "nav_dead_1.2s_gyro_only.csv"},

        // 第二组：陀螺+加计全组合 (压制舒拉震荡)
        //{"Group_2.0s_GyroAcc",  2.0, 0.0, true,  "nav_dead_2.0s_gyro_acc.csv"}, 
        {"Group_1.6s_GyroAcc",  1.6, 0.4, true,  "nav_dead_1.6s_gyro_acc.csv"}, 
        //{"Group_1.2s_GyroAcc",  1.2, 0.8, true,  "nav_dead_1.2s_gyro_acc.csv"}
    };

    cout << "\n>>> Starting Batch Processing (" << experiments.size() << " scenarios)..." << endl;
    
    for (const auto& exp : experiments) {
        // 调用封装好的函数
        RunSingleExperiment(exp, ctx);
    }

    cout << "\nAll Scenarios Completed Successfully!" << endl;
    return 0;
}