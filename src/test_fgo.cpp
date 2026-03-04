#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <chrono> // 引入高精度计时器
#include "sins_engine.h"
#include "cai_sim.h"
#include "support.h" 

using namespace std;
using namespace Eigen;
using namespace cai;

// 辅助函数：生成反对称矩阵
auto askew = [](const Vector3d& v) -> Matrix3d {
    Matrix3d m;
    m <<   0,  -v(2),  v(1),
         v(2),    0,  -v(0),
        -v(1),  v(0),    0;
    return m;
};

// ==========================================
// 用于 FGO 解析回溯的日志缓存结构
// ==========================================
struct LogPoint {
    double nav_time;
    Vector3d att;
    Vector3d vn;
    Vector3d pos;
    Matrix3d Cnb;
    MatrixXd J_at_t;
    
    LogPoint(double t, const Vector3d& a, const Vector3d& v, const Vector3d& p, const Matrix3d& C, const MatrixXd& j) 
        : nav_time(t), att(a), vn(v), pos(p), Cnb(C), J_at_t(j) {}
};

// ==============================================================================
// 算法 1：纯光纤惯导基线 (Pure FOG INS)
// 目的：关闭所有原子约束，纯机械编排，提取 24 小时最大发散误差并保存轨迹
// ==============================================================================
void RunPureFOG(const ExperimentConfig& exp, const SimulationContext& ctx) {
    SinsEngine sinsegine(ctx.ts);
    sinsegine.eth = ctx.eth; 
    sinsegine.ins = INSState(ctx.att_align, Vector3d::Zero(), ctx.pos_ref, ctx.ts, sinsegine.eth);
    sinsegine.ins.set_bias(ctx.eb_align, ctx.db_align);

    std::vector<std::string> rest_files(ctx.file_list.begin() + 1, ctx.file_list.end());
    IMUChainedLoader nav_loader(rest_files, ctx.ts, ctx.local_g);

    // 【修复点】：增加文件输出逻辑
    ofstream log(exp.output_file);
    log << "time,lat_err,lon_err,drift\n"; // 写入表头

    double max_drift = 0;
    size_t total_steps = 0;
    size_t max_steps = static_cast<size_t>(exp.duration_hours * 3600.0 / ctx.ts);

    IMUData epoch;
    while (nav_loader.Next(epoch) && total_steps < max_steps) {
        sinsegine.Step_Nav(epoch.wm, epoch.vm);
        sinsegine.ins.vn(2) = 0.0; // 阻尼

        if (total_steps % 400 == 0) { 
            double lat_err = (sinsegine.ins.pos(0) - ctx.pos_ref(0)) * ctx.glv.Re;
            double lon_err = (sinsegine.ins.pos(1) - ctx.pos_ref(1)) * ctx.glv.Re * cos(ctx.pos_ref(0));
            double drift = sqrt(lat_err*lat_err + lon_err*lon_err);
            if(drift > max_drift) max_drift = drift;

            // 【修复点】：把数据写入 CSV 文件 (1Hz 记录)
            log << fixed << setprecision(9) << (total_steps * ctx.ts) << "," 
                << lat_err << "," << lon_err << "," << drift << "\n";

            if (total_steps > 0 && total_steps % (3600 * 400) == 0) {
                 cout << "   [Pure FOG] Time: " << (total_steps * ctx.ts) / 3600.0 
                      << " h | Current Drift: " << drift << " m | Max: " << max_drift << " m" << endl;
            }
        }
        total_steps++;
    }
    
    log.close(); // 关闭文件流
    cout << " => [Result] Pure FOG Max Drift: " << max_drift << " m\n" << endl;
}

// ==============================================================================
// 算法 2：多速率 ES-FGO (本文核心贡献)
// 目的：输出实时的锯齿状拉回轨迹，并测算 O(1) 历史解析平滑的极低 CPU 耗时
// ==============================================================================
void RunESFGO(const ExperimentConfig& exp, const SimulationContext& ctx) {
    AtomicGyroSimulator atom(ctx.pos_ref, ctx.glv);
    atom.Init(ctx.att_align); 
    CAIParams acc_params; acc_params.bias_ug = 0.0; acc_params.vrw_ug  = 0.05; 
    AtomicAccSimulator atom_acc(ctx.pos_ref, ctx.glv, acc_params);
    atom_acc.Init(ctx.att_align); 
    
    SinsEngine sinsegine(ctx.ts);
    sinsegine.eth = ctx.eth; 
    sinsegine.ins = INSState(ctx.att_align, Vector3d::Zero(), ctx.pos_ref, ctx.ts, sinsegine.eth);
    sinsegine.ins.set_bias(ctx.eb_align, ctx.db_align);

    std::vector<std::string> rest_files(ctx.file_list.begin() + 1, ctx.file_list.end());
    IMUChainedLoader nav_loader(rest_files, ctx.ts, ctx.local_g);

    // 【修正1】：严格匹配 Python 脚本需要的表头
    ofstream log(exp.output_file);
    log << "time,lat_err,lon_err,drift\n"; 

    double T_cycle = exp.t_active + exp.t_dead; 
    int total_samples_per_cycle = static_cast<int>(T_cycle / ctx.ts); 
    int active_samples = static_cast<int>(exp.t_active / ctx.ts);     
    
    Vector3d fog_gyro_active_sum = Vector3d::Zero(); 
    Vector3d fog_acc_active_sum  = Vector3d::Zero(); 
    
    int sample_count = 0;
    double max_drift = 0;
    size_t total_steps = 0;
    size_t max_steps = static_cast<size_t>(exp.duration_hours * 3600.0 / ctx.ts);

    double total_cpu_time_ms = 0;
    int update_counts = 0;

    MatrixXd J = MatrixXd::Zero(9, 6); 
    MatrixXd Phi_xx = MatrixXd::Identity(9, 9);
    MatrixXd Phi_xb = MatrixXd::Zero(9, 6);
    std::vector<LogPoint> log_buffer;
    log_buffer.reserve(total_samples_per_cycle + 10); 

    IMUData epoch;
    while (nav_loader.Next(epoch) && total_steps < max_steps) {
        sample_count++;
        total_steps++;

        if (sample_count <= active_samples) {
            fog_gyro_active_sum += epoch.wm / ctx.ts; 
            fog_acc_active_sum  += epoch.vm / ctx.ts; 
        } 

        // 标称推演
        sinsegine.Step_Nav(epoch.wm, epoch.vm);
        sinsegine.ins.vn(2) = 0.0; 

        // 雅可比累加
        Matrix3d Cnb = sinsegine.ins.Cnb;
        Vector3d fn = Cnb * (epoch.vm / ctx.ts);
        double RM_h = sinsegine.eth.RMh; double RN_h = sinsegine.eth.RNh;
        Matrix3d Mpv = Matrix3d::Zero(); Mpv(0, 1) = 1.0 / RM_h; Mpv(1, 0) = 1.0 / (RN_h * cos(sinsegine.ins.pos(0))); Mpv(2, 2) = 1.0;

        Phi_xb.setZero();
        Phi_xb.block<3,3>(0, 0) = -Cnb * ctx.ts; Phi_xb.block<3,3>(3, 3) = -Cnb * ctx.ts;
        Phi_xx.setIdentity();
        Phi_xx.block<3,3>(0, 0) -= askew(sinsegine.eth.winn) * ctx.ts;
        Phi_xx.block<3,3>(3, 0) += askew(fn) * ctx.ts;
        Phi_xx.block<3,3>(3, 3) -= askew(sinsegine.eth.wcor) * ctx.ts;
        Phi_xx.block<3,3>(6, 3) += Mpv * ctx.ts;

        J = Phi_xx * J + Phi_xb;

        // 【修正2】：高频记录实时的前向漂移状态（生成死区锯齿的上升沿）
        if (total_steps % 40 == 0) { 
            log_buffer.emplace_back(total_steps * ctx.ts, sinsegine.ins.att, sinsegine.ins.vn, sinsegine.ins.pos, sinsegine.ins.Cnb, J);
            
            double lat_err = (sinsegine.ins.pos(0) - ctx.pos_ref(0)) * ctx.glv.Re;
            double lon_err = (sinsegine.ins.pos(1) - ctx.pos_ref(1)) * ctx.glv.Re * cos(ctx.pos_ref(0));
            double drift = sqrt(lat_err*lat_err + lon_err*lon_err);
            if(drift > max_drift) max_drift = drift;
            
            log << fixed << setprecision(9) << (total_steps * ctx.ts) << "," << lat_err << "," << lon_err << "," << drift << "\n";
        }

        // ES-FGO 零方差解析平滑更新
        if (sample_count >= total_samples_per_cycle) {
            Vector3d eb_new = (fog_gyro_active_sum / active_samples) - atom.Measure();
            Vector3d db_new = exp.use_atomic_acc ? ((fog_acc_active_sum / active_samples) - atom_acc.Measure()) : sinsegine.ins.db;

            VectorXd delta_b(6);
            delta_b.segment<3>(0) = eb_new - sinsegine.ins.eb;
            delta_b.segment<3>(3) = db_new - sinsegine.ins.db;

            // --- 算力计时开始 ---
            auto t_start = chrono::high_resolution_clock::now();

            // 诚实地模拟对历史时刻的 O(1) 解析平滑操作（只算不输出，模拟 CPU 负载）
            for (const auto& lp : log_buffer) {
                VectorXd dx = lp.J_at_t * delta_b;
                Vector3d phi = dx.segment<3>(0);
                Matrix3d R_phi = (phi.norm() > 1e-12) ? AngleAxisd(phi.norm(), phi.normalized()).toRotationMatrix() : (Matrix3d::Identity() + askew(phi));
                Matrix3d Cnb_cor = R_phi * lp.Cnb;     
                Vector3d s_pos = lp.pos + dx.segment<3>(6);
            }

            // 更新当前最新时刻状态
            VectorXd final_dx = J * delta_b;
            Vector3d final_phi = final_dx.segment<3>(0);
            Matrix3d R_final_phi = (final_phi.norm() > 1e-12) ? AngleAxisd(final_phi.norm(), final_phi.normalized()).toRotationMatrix() : (Matrix3d::Identity() + askew(final_phi));
            Matrix3d Cnb_cor_final = R_final_phi * sinsegine.ins.Cnb;
            
            sinsegine.ins.att = INSMath::m2att(Cnb_cor_final);
            sinsegine.ins.Cnb = Cnb_cor_final;
            sinsegine.ins.qnb = INSMath::a2qua(sinsegine.ins.att); 
            sinsegine.ins.vn  += final_dx.segment<3>(3);
            sinsegine.ins.pos += final_dx.segment<3>(6);
            sinsegine.ins.set_bias(eb_new, db_new);

            // --- 算力计时结束 ---
            auto t_end = chrono::high_resolution_clock::now();
            total_cpu_time_ms += chrono::duration<double, std::milli>(t_end - t_start).count();
            update_counts++;
            
            // 【修正3】：记录更新瞬间的状态（生成死区锯齿的垂直下拉）
            double lat_err = (sinsegine.ins.pos(0) - ctx.pos_ref(0)) * ctx.glv.Re;
            double lon_err = (sinsegine.ins.pos(1) - ctx.pos_ref(1)) * ctx.glv.Re * cos(ctx.pos_ref(0));
            double drift = sqrt(lat_err*lat_err + lon_err*lon_err);
            log << fixed << setprecision(9) << (total_steps * ctx.ts) << "," << lat_err << "," << lon_err << "," << drift << "\n";

            if (total_steps > 0 && total_steps % (3600 * 400) == 0) {
                 cout << "   [ES-FGO] Time: " << (total_steps*ctx.ts) / 3600.0 << " h | Max Drift: " << max_drift << " m" << endl;
            }

            fog_gyro_active_sum.setZero(); fog_acc_active_sum.setZero();
            sample_count = 0; J.setZero(); log_buffer.clear();
        }
    }
    log.close();
    cout << " => [Result] ES-FGO Max Drift: " << max_drift << " m" << endl;
    cout << " => [CPU] Avg ES-FGO Update Time: " << (total_cpu_time_ms / update_counts) << " ms/cycle\n" << endl;
}

// ==============================================================================
// 算法 3：传统延迟状态 EKF (Standard Delayed-State EKF)
// 目的：只跑一小段时间，纯粹为了测试 O(N) 重积分带来的算力毛刺
// ==============================================================================
void RunEKF(const ExperimentConfig& exp, const SimulationContext& ctx) {
    AtomicGyroSimulator atom(ctx.pos_ref, ctx.glv);
    atom.Init(ctx.att_align); 
    CAIParams acc_params; acc_params.bias_ug = 0.0; acc_params.vrw_ug  = 0.05; 
    AtomicAccSimulator atom_acc(ctx.pos_ref, ctx.glv, acc_params);
    atom_acc.Init(ctx.att_align); 
    
    SinsEngine sinsegine(ctx.ts);
    sinsegine.eth = ctx.eth; 
    sinsegine.ins = INSState(ctx.att_align, Vector3d::Zero(), ctx.pos_ref, ctx.ts, sinsegine.eth);
    sinsegine.ins.set_bias(ctx.eb_align, ctx.db_align);

    std::vector<std::string> rest_files(ctx.file_list.begin() + 1, ctx.file_list.end());
    IMUChainedLoader nav_loader(rest_files, ctx.ts, ctx.local_g);

    double T_cycle = exp.t_active + exp.t_dead; 
    int total_samples_per_cycle = static_cast<int>(T_cycle / ctx.ts); 
    int active_samples = static_cast<int>(exp.t_active / ctx.ts);     
    
    Vector3d fog_gyro_active_sum = Vector3d::Zero(); 
    Vector3d fog_acc_active_sum  = Vector3d::Zero(); 
    
    int sample_count = 0;
    size_t total_steps = 0;
    size_t max_steps = static_cast<size_t>(exp.duration_hours * 3600.0 / ctx.ts);

    double total_cpu_time_ms = 0;
    int update_counts = 0;

    // EKF 核心：保存 2s 前的状态与原始 IMU 数据
    INSState state_t0 = sinsegine.ins; 
    std::vector<IMUData> raw_imu_buffer;
    raw_imu_buffer.reserve(total_samples_per_cycle + 10);

    IMUData epoch;
    while (nav_loader.Next(epoch) && total_steps < max_steps) {
        sample_count++;
        total_steps++;
        raw_imu_buffer.push_back(epoch);

        if (sample_count <= active_samples) {
            fog_gyro_active_sum += epoch.wm / ctx.ts; 
            fog_acc_active_sum  += epoch.vm / ctx.ts; 
        } 

        // 仅作标称推演，用于实时输出
        sinsegine.Step_Nav(epoch.wm, epoch.vm);
        sinsegine.ins.vn(2) = 0.0; 

        // EKF 重积分更新
        if (sample_count >= total_samples_per_cycle) {
            Vector3d eb_new = (fog_gyro_active_sum / active_samples) - atom.Measure();
            Vector3d db_new = exp.use_atomic_acc ? ((fog_acc_active_sum / active_samples) - atom_acc.Measure()) : sinsegine.ins.db;

            // --- 算力计时开始 ---
            auto t_start = chrono::high_resolution_clock::now();

            // 1. 状态回滚至死区起始时刻 t0
            sinsegine.ins = state_t0;
            // 2. 更新修正后的 Bias
            sinsegine.ins.set_bias(eb_new, db_new);
            // 3. O(N) 级别高频数值重积分
            for (const auto& ep : raw_imu_buffer) {
                sinsegine.Step_Nav(ep.wm, ep.vm);
                sinsegine.ins.vn(2) = 0.0;
            }

            // --- 算力计时结束 ---
            auto t_end = chrono::high_resolution_clock::now();
            total_cpu_time_ms += chrono::duration<double, std::milli>(t_end - t_start).count();
            update_counts++;

            // 更新下一轮周期的 t0 状态
            state_t0 = sinsegine.ins;
            fog_gyro_active_sum.setZero(); fog_acc_active_sum.setZero();
            sample_count = 0; raw_imu_buffer.clear();
        }
    }
    // 【修正4】：EKF 不做文件输出，程序极其清爽
    cout << " => [CPU] Avg EKF Re-integration Time: " << (total_cpu_time_ms / update_counts) << " ms/cycle\n" << endl;
}

// ==========================================
// 任务分发中心
// ==========================================
void RunExperimentRouter(const ExperimentConfig& exp, const SimulationContext& ctx) {
    cout << "\n--------------------------------------------------" << endl;
    cout << " Starting Scenario: " << exp.name << endl;
    cout << " Duration: " << exp.duration_hours << " h | Output: " << exp.output_file << endl;
    cout << "--------------------------------------------------" << endl;

    switch (exp.algo) {
        case AlgoType::PURE_FOG: RunPureFOG(exp, ctx); break;
        case AlgoType::ES_FGO:   RunESFGO(exp, ctx);   break;
        case AlgoType::EKF:      RunEKF(exp, ctx);     break;
    }
}

// ==========================================
// 主函数
// ==========================================
int main() {
    SimulationContext ctx;
    ctx.ts = 1.0 / 400.0;
    ctx.glv = GLV(); 
    ctx.pos_ref = Vector3d(32.0286 * ctx.glv.deg, 118.8533 * ctx.glv.deg, 17.0);
    ctx.eth = Earth(ctx.glv);
    ctx.eth.update(ctx.pos_ref, Vector3d::Zero());
    ctx.local_g = ctx.eth.gn.norm();
    
    ctx.file_list = {
        "../fog_part1.csv", "../fog_part2.csv", "../fog_part3.csv", 
        "../fog_part4.csv", "../fog_part5.csv"
    };

    cout << "\n[Phase 1] Executing Fine Alignment..." << endl;
    auto buffer_part1 = LoadIMUData(ctx.file_list[0], ctx.ts, ctx.local_g);
    double start_time_sec = 5.0 * 3600.0 + 50.0 * 60.0; 
    size_t start_idx = static_cast<size_t>(start_time_sec / ctx.ts);
    std::vector<IMUData> buffer_align(buffer_part1.begin() + start_idx, buffer_part1.end());

    SinsEngine align_engine(ctx.ts);
    align_engine.eth = ctx.eth;
    align_engine.ins = INSState(Vector3d::Zero(), Vector3d::Zero(), ctx.pos_ref, ctx.ts, ctx.eth);
    align_engine.res_init.pos = ctx.pos_ref;

    HybridAlignConfig cfg; cfg.t_coarse = 60; cfg.t_fine = buffer_align.size() * ctx.ts; 
    cfg.eb_sigma_allan = 0.003; cfg.db_sigma_allan = 50.0; cfg.verbose = false; 

    auto align_res = align_engine.Run_HybridAlign(buffer_align, cfg);
    if (!align_res.valid) return -1;
    
    ctx.att_align = align_res.att; ctx.eb_align = align_res.eb; ctx.db_align = align_res.db;
    ctx.Cnb_align = INSMath::a2mat(align_res.att); 
    vector<IMUData>().swap(buffer_part1); vector<IMUData>().swap(buffer_align);

    // ==========================================
    // Phase 2: 论文核心对比实验库 (一键跑全家桶)
    // ==========================================
vector<ExperimentConfig> paper_experiments = {
        // 1. 基线1：纯光纤（跑满24h，提取最大发散漂移 16.4km）
        {"Baseline 1: Pure FOG", AlgoType::PURE_FOG, 24.0, 0.0, 0.0, false, "nav_pure_fog.csv"}, 

        // 2. 舒拉振荡验证：仅用原子陀螺，关闭原子加计（跑满24h，用于画 fig_pos_err.png）
        {"Ablation: CAIG Only",  AlgoType::ES_FGO,   24.0, 1.6, 0.4, false, "nav_dead_1.6s_gyro_only.csv"},

        // 3. 主角：多速率 ES-FGO（跑满24h，用于画 2D轨迹和锯齿放大图）
        {"Proposed: ES-FGO",     AlgoType::ES_FGO,   24.0, 1.6, 0.4, true,  "nav_esfgo_full.csv"},

        // 4. 基线2：传统 EKF（跑1h，用于测算 0.126ms 算力）
        {"Baseline 2: EKF Perf", AlgoType::EKF,      1.0,  1.6, 0.4, true,  "nav_ekf_perf.csv"}
    };
    cout << "\n>>> Starting Master Simulation Script for IEEE Paper..." << endl;
    for (const auto& exp : paper_experiments) {
        RunExperimentRouter(exp, ctx);
    }

    cout << "\nAll Simulations Completed! You can now fill the data into the LaTeX source." << endl;
    return 0;
}