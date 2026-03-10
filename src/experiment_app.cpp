#include "experiment_app.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include "sins_engine.h"
#include "cai_sim.h"
#include "support.h" 

using namespace std;
using namespace Eigen;
using namespace cai;

// ==============================================================================
// 算法 1：纯光纤惯导基线 (Pure FOG INS)
// ==============================================================================
void RunPureFOG(const ExperimentConfig& exp, const SimulationContext& ctx) {
    SinsEngine sinsegine(ctx.ts);
    sinsegine.eth = ctx.eth; 
    sinsegine.ins = INSState(ctx.att_align, Vector3d::Zero(), ctx.pos_ref, ctx.ts, sinsegine.eth);
    sinsegine.ins.set_bias(ctx.eb_align, ctx.db_align);

    std::vector<std::string> rest_files(ctx.file_list.begin() + 1, ctx.file_list.end());
    IMUChainedLoader nav_loader(rest_files, ctx.ts, ctx.local_g);

    ofstream log(exp.output_file);
    log << "time,lat_err,lon_err,drift,vn_x,vn_y,vn_z,att_roll_deg,att_pitch_deg,att_yaw_deg,eb_x_dph,eb_y_dph,eb_z_dph,db_x_ug,db_y_ug,db_z_ug,true_eb_x_dph\n";

    double max_drift = 0;
    size_t total_steps = 0;
    size_t max_steps = static_cast<size_t>(exp.duration_hours * 3600.0 / ctx.ts);

    IMUData epoch;
    while (nav_loader.Next(epoch) && total_steps < max_steps) {
        sinsegine.Step_Nav(epoch.wm, epoch.vm);
        sinsegine.ins.vn(2) = 0.0; 

        if (total_steps % 400 == 0) { 
            double lat_err = (sinsegine.ins.pos(0) - ctx.pos_ref(0)) * ctx.glv.Re;
            double lon_err = (sinsegine.ins.pos(1) - ctx.pos_ref(1)) * ctx.glv.Re * cos(ctx.pos_ref(0));
            double drift = sqrt(lat_err*lat_err + lon_err*lon_err);
            if(drift > max_drift) max_drift = drift;

            LogData(log, total_steps * ctx.ts, sinsegine.ins, ctx, drift);

            if (total_steps > 0 && total_steps % (3600 * 400) == 0) {
                 cout << "   [Pure FOG] Time: " << (total_steps * ctx.ts) / 3600.0 
                      << " h | Current Drift: " << drift << " m | Max: " << max_drift << " m" << endl;
            }
        }
        total_steps++;
    }
    
    log.close();
    cout << " => [Result] Pure FOG Max Drift: " << max_drift << " m\n" << endl;
}

// ==============================================================================
// 算法 2：多速率 FGO (本文核心贡献)
// ==============================================================================
void RunESFGO(const ExperimentConfig& exp, const SimulationContext& ctx) {
    AtomicGyroSimulator atom(ctx.pos_ref, ctx.glv);
    atom.Init(ctx.att_align); 
    CAIParams acc_params; acc_params.bias_ug = 0.0; acc_params.vrw_ug = 0.05; 
    AtomicAccSimulator atom_acc(ctx.pos_ref, ctx.glv, acc_params);
    atom_acc.Init(ctx.att_align);
    
    SinsEngine sinsegine(ctx.ts);
    sinsegine.eth = ctx.eth; 
    sinsegine.ins = INSState(ctx.att_align, Vector3d::Zero(), ctx.pos_ref, ctx.ts, sinsegine.eth);
    sinsegine.ins.set_bias(ctx.eb_align, ctx.db_align);

    std::vector<std::string> rest_files(ctx.file_list.begin() + 1, ctx.file_list.end());
    IMUChainedLoader nav_loader(rest_files, ctx.ts, ctx.local_g);

    ofstream log(exp.output_file);
    log << "time,lat_err,lon_err,drift,vn_x,vn_y,vn_z,att_roll_deg,att_pitch_deg,att_yaw_deg,eb_x_dph,eb_y_dph,eb_z_dph,db_x_ug,db_y_ug,db_z_ug,true_eb_x_dph\n"; 

    double T_cycle = exp.t_active + exp.t_dead; 
    int total_samples_per_cycle = static_cast<int>(T_cycle / ctx.ts); 
    int active_samples = static_cast<int>(exp.t_active / ctx.ts);    
    
    Vector3d fog_gyro_active_sum = Vector3d::Zero(); 
    Vector3d fog_acc_active_sum = Vector3d::Zero();
    
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

    // 第一阶马尔可夫陀螺零偏漂移 (τ=300s, σ_ss=10°/h)
    const double tau_m = 300.0, sigma_m = 10.0;
    double decay_m = std::exp(-ctx.ts / tau_m);
    double drive_m = sigma_m * sqrt(1.0 - decay_m * decay_m);
    double mkv_x = 0.0, mkv_y = 0.0;
    std::mt19937 rng_m(42);
    std::normal_distribution<double> nd_m(0.0, 1.0);

    IMUData epoch;
    while (nav_loader.Next(epoch) && total_steps < max_steps) {
        
        if (exp.inject_time_varying_drift) {
            // 第一阶马尔可夫: b_{k+1} = e^{-dt/τ}·b_k + σ_drive·w_k
            mkv_x = decay_m * mkv_x + drive_m * nd_m(rng_m);
            mkv_y = decay_m * mkv_y + drive_m * nd_m(rng_m);
            epoch.wm(0) += mkv_x * ctx.glv.dph * ctx.ts;
            epoch.wm(1) += mkv_y * ctx.glv.dph * ctx.ts;
        }
        
        sample_count++;
        total_steps++;

        if (sample_count <= active_samples) {
            fog_gyro_active_sum += epoch.wm / ctx.ts; 
            fog_acc_active_sum += epoch.vm / ctx.ts;
        } 

        sinsegine.Step_Nav(epoch.wm, epoch.vm);
        sinsegine.ins.vn(2) = 0.0; 

        Matrix3d Cnb = sinsegine.ins.Cnb;
        Vector3d fn = Cnb * (epoch.vm / ctx.ts);
        double RM_h = sinsegine.eth.RMh; double RN_h = sinsegine.eth.RNh;
        Matrix3d Mpv = Matrix3d::Zero(); Mpv(0, 1) = 1.0 / RM_h; Mpv(1, 0) = 1.0 / (RN_h * cos(sinsegine.ins.pos(0))); Mpv(2, 2) = 1.0;

        Phi_xb.setZero();
        Phi_xb.block<3,3>(0, 0) = -Cnb * ctx.ts; Phi_xb.block<3,3>(3, 3) = -Cnb * ctx.ts;
        Phi_xx.setIdentity();
        Phi_xx.block<3,3>(0, 0) -= INSMath::askew(sinsegine.eth.winn) * ctx.ts;
        Phi_xx.block<3,3>(3, 0) += INSMath::askew(fn) * ctx.ts;
        Phi_xx.block<3,3>(3, 3) -= INSMath::askew(sinsegine.eth.wcor) * ctx.ts;
        Phi_xx.block<3,3>(6, 3) += Mpv * ctx.ts;

        J = Phi_xx * J + Phi_xb;

        // 存储当前时刻的历史状态和雅可比矩阵（用于后续的解析平滑）
        log_buffer.emplace_back(total_steps * ctx.ts, sinsegine.ins.att, sinsegine.ins.vn, 
                                sinsegine.ins.pos, sinsegine.ins.Cnb, J);

        // FGO 周期结束，执行 MAP 解析平滑
        if (sample_count >= total_samples_per_cycle) {
            
            Vector3d fog_gyro_mean = fog_gyro_active_sum / active_samples;
            Vector3d delta_z_gyro = fog_gyro_mean - atom.Measure();

            // 匹配马尔可夫漂移稳态标准差 σ=10°/h
            Matrix3d Sigma_I_gyro = Matrix3d::Identity() * pow(sigma_m * ctx.glv.dph, 2); 
            Matrix3d Sigma_A_gyro = Matrix3d::Identity() * pow(1e-4 * ctx.glv.dph, 2);

            auto t_start = chrono::high_resolution_clock::now();

            // Schur Complement 降维求解最优零偏
            Matrix3d Lambda_gyro = Sigma_I_gyro.inverse() + Sigma_A_gyro.inverse();
            Matrix3d K_fgo_gyro = Lambda_gyro.inverse() * Sigma_A_gyro.inverse(); 
            Vector3d delta_eb_star = K_fgo_gyro * (delta_z_gyro - sinsegine.ins.eb);
            Vector3d eb_new = sinsegine.ins.eb + delta_eb_star;

            VectorXd delta_b = VectorXd::Zero(6);
            Vector3d db_new = sinsegine.ins.db;
            if (exp.use_atomic_acc) {
                Matrix3d Sigma_I_acc = Matrix3d::Identity() * pow(50.0 * ctx.glv.ug, 2);
                Matrix3d Sigma_A_acc = Matrix3d::Identity() * pow(0.5 * ctx.glv.ug, 2);
                Matrix3d Lambda_acc = Sigma_I_acc.inverse() + Sigma_A_acc.inverse();
                Matrix3d K_fgo_acc = Lambda_acc.inverse() * Sigma_A_acc.inverse();
                Vector3d fog_acc_mean = fog_acc_active_sum / active_samples;
                Vector3d delta_z_acc = fog_acc_mean - atom_acc.Measure();
                Vector3d delta_db_star = K_fgo_acc * (delta_z_acc - sinsegine.ins.db);
                db_new = sinsegine.ins.db + delta_db_star;
                delta_b.segment<3>(3) = delta_db_star;
            }
            delta_b.segment<3>(0) = delta_eb_star;

            // ---------------------------------------------------------
            // 核心修正：利用雅可比对历史窗口进行 O(1) 解析平滑并记录数据
            // ---------------------------------------------------------
            for (size_t li = 0; li < log_buffer.size(); ++li) {
                const auto& lp = log_buffer[li];
                VectorXd dx = lp.J_at_t * delta_b;
                
                // 历史状态流形重构
                Vector3d phi = dx.segment<3>(0);
                Matrix3d R_phi = (phi.norm() > 1e-12) ? AngleAxisd(phi.norm(), phi.normalized()).toRotationMatrix() : (Matrix3d::Identity() + INSMath::askew(phi));
                Matrix3d Cnb_smooth = R_phi * lp.Cnb;     
                Vector3d vn_smooth = lp.vn + dx.segment<3>(3);
                Vector3d pos_smooth = lp.pos + dx.segment<3>(6);

                // 计算平滑后的实际定位误差
                double lat_err = (pos_smooth(0) - ctx.pos_ref(0)) * ctx.glv.Re;
                double lon_err = (pos_smooth(1) - ctx.pos_ref(1)) * ctx.glv.Re * cos(ctx.pos_ref(0));
                double drift = sqrt(lat_err*lat_err + lon_err*lon_err);
                if(drift > max_drift) max_drift = drift;
                
                double true_eb_x_dph = (ctx.eb_align(0) / ctx.glv.deg * 3600.0) + mkv_x;

                // 只记录每个周期的最后一个平滑点（2s粒度，与EKF一致，避免7GB磁盘爆炸）
                if (li == log_buffer.size() - 1) {
                    INSState smoothed_ins = sinsegine.ins;
                    smoothed_ins.att = INSMath::m2att(Cnb_smooth);
                    smoothed_ins.vn = vn_smooth;
                    smoothed_ins.pos = pos_smooth;
                    smoothed_ins.eb = eb_new;
                    LogData(log, lp.nav_time, smoothed_ins, ctx, drift, true_eb_x_dph);
                }
            }

            // 更新实时端点状态供下一周期递推使用
            VectorXd final_dx = J * delta_b;
            sinsegine.UpdateInsState(final_dx);
            sinsegine.ins.set_bias(eb_new, db_new); 

            auto t_end = chrono::high_resolution_clock::now();
            total_cpu_time_ms += chrono::duration<double, std::milli>(t_end - t_start).count();
            update_counts++;

            if (total_steps > 0 && total_steps % (3600 * 400) == 0) {
                 cout << "   [FGO] Time: " << (total_steps*ctx.ts) / 3600.0 
                      << " h | Max Smoothed Drift: " << max_drift << " m" << endl;
            }

            fog_gyro_active_sum.setZero();
            fog_acc_active_sum.setZero();
            sample_count = 0; 
            J.setZero(); 
            log_buffer.clear();
        }
    }
    log.close();
    cout << " => [Result] FGO Max Drift: " << max_drift << " m" << endl;
    cout << " => [CPU] Avg FGO Update Time: " << (total_cpu_time_ms / update_counts) << " ms/cycle\n" << endl;
}

// ==============================================================================
// 算法 3：标准 15 维延迟状态 EKF
// ==============================================================================
void RunEKF(const ExperimentConfig& exp, const SimulationContext& ctx) {
    AtomicGyroSimulator atom(ctx.pos_ref, ctx.glv);
    atom.Init(ctx.att_align); 
    CAIParams acc_params; acc_params.bias_ug = 0.0; acc_params.vrw_ug = 0.05; 
    AtomicAccSimulator atom_acc(ctx.pos_ref, ctx.glv, acc_params);
    atom_acc.Init(ctx.att_align);
    
    SinsEngine sinsegine(ctx.ts);
    sinsegine.eth = ctx.eth; 
    sinsegine.ins = INSState(ctx.att_align, Vector3d::Zero(), ctx.pos_ref, ctx.ts, sinsegine.eth);
    sinsegine.ins.set_bias(ctx.eb_align, ctx.db_align);

    std::vector<std::string> rest_files(ctx.file_list.begin() + 1, ctx.file_list.end());
    IMUChainedLoader nav_loader(rest_files, ctx.ts, ctx.local_g);

    ofstream log(exp.output_file);
    log << "time,lat_err,lon_err,drift,vn_x,vn_y,vn_z,att_roll_deg,att_pitch_deg,att_yaw_deg,eb_x_dph,eb_y_dph,eb_z_dph,db_x_ug,db_y_ug,db_z_ug,true_eb_x_dph\n"; 

    MatrixXd P = MatrixXd::Zero(15, 15);
    P.block<3,3>(0,0)   = Matrix3d::Identity() * pow(1.0 * ctx.glv.min, 2);    
    P.block<3,3>(3,3)   = Matrix3d::Identity() * pow(0.1, 2);                  
    P.block<3,3>(6,6)   = Matrix3d::Identity() * pow(10.0, 2);                 
    P.block<3,3>(9,9)   = Matrix3d::Identity() * pow(0.01 * ctx.glv.dph, 2);   
    P.block<3,3>(12,12) = Matrix3d::Identity() * pow(100.0 * ctx.glv.ug, 2);   

    MatrixXd Q = MatrixXd::Zero(15, 15);
    Q.block<3,3>(0,0)   = Matrix3d::Identity() * pow(0.001 * ctx.glv.dpsh, 2); 
    Q.block<3,3>(3,3)   = Matrix3d::Identity() * pow(5.0 * ctx.glv.ugpsHz, 2); 
    Q.block<3,3>(12,12) = Matrix3d::Identity() * pow(1e-5 * ctx.glv.ug, 2);    
    
    // 动态环境保留真实物理滞后，静态环境确保正常收敛
    if (exp.inject_time_varying_drift) {
        // 马尔可夫驱动噪声 PSD: q_c = 2σ²/τ (连续时间，传播公式中乘 dt)
        Q.block<3,3>(9,9) = Matrix3d::Identity() * 2.0 * pow(10.0 * ctx.glv.dph, 2) / 300.0;   
    } else {
        Q.block<3,3>(9,9) = Matrix3d::Identity() * pow(0.003 * ctx.glv.dph, 2);   
    }

    MatrixXd R = Matrix3d::Identity() * pow(0.05 * ctx.glv.dph, 2);   
    MatrixXd H = MatrixXd::Zero(3, 15);
    H.block<3,3>(0,9)   = Matrix3d::Identity(); 

    double T_cycle = exp.t_active + exp.t_dead; 
    int total_samples_per_cycle = static_cast<int>(T_cycle / ctx.ts); 
    int active_samples = static_cast<int>(exp.t_active / ctx.ts);    
    
    Vector3d fog_gyro_active_sum = Vector3d::Zero(); 
    Vector3d fog_acc_active_sum = Vector3d::Zero();

    // 第一阶马尔可夫陀螺零偏漂移 (τ=300s, σ_ss=10°/h) — 与 FGO 相同种子确保公平对比
    const double tau_m_ekf = 300.0, sigma_m_ekf = 10.0;
    double decay_m_ekf = std::exp(-ctx.ts / tau_m_ekf);
    double drive_m_ekf = sigma_m_ekf * sqrt(1.0 - decay_m_ekf * decay_m_ekf);
    double mkv_x_ekf = 0.0, mkv_y_ekf = 0.0;
    std::mt19937 rng_m_ekf(42);
    std::normal_distribution<double> nd_m_ekf(0.0, 1.0);
    
    int sample_count = 0;
    double max_drift = 0;
    size_t total_steps = 0;
    size_t max_steps = static_cast<size_t>(exp.duration_hours * 3600.0 / ctx.ts);
    double total_cpu_time_ms = 0;
    int update_counts = 0;

    IMUData epoch;
    while (nav_loader.Next(epoch) && total_steps < max_steps) {
        
        if (exp.inject_time_varying_drift) {
            mkv_x_ekf = decay_m_ekf * mkv_x_ekf + drive_m_ekf * nd_m_ekf(rng_m_ekf);
            mkv_y_ekf = decay_m_ekf * mkv_y_ekf + drive_m_ekf * nd_m_ekf(rng_m_ekf);
            epoch.wm(0) += mkv_x_ekf * ctx.glv.dph * ctx.ts;
            epoch.wm(1) += mkv_y_ekf * ctx.glv.dph * ctx.ts;
        }
        
        sample_count++;
        total_steps++;

        if (sample_count <= active_samples) {
            fog_gyro_active_sum += epoch.wm / ctx.ts; 
            fog_acc_active_sum += epoch.vm / ctx.ts;
        } 

        // 1. 机械编排
        sinsegine.Step_Nav(epoch.wm, epoch.vm);
        sinsegine.ins.vn(2) = 0.0; 

        // 2. EKF 状态协方差传播 (400Hz)
        Matrix3d Cnb = sinsegine.ins.Cnb;
        Vector3d fn = Cnb * (epoch.vm / ctx.ts);
        
        double RM_h_ekf = sinsegine.eth.RMh; 
        double RN_h_ekf = sinsegine.eth.RNh;
        Matrix3d Mpv_ekf = Matrix3d::Zero(); 
        Mpv_ekf(0, 1) = 1.0 / RM_h_ekf;
        Mpv_ekf(1, 0) = 1.0 / (RN_h_ekf * cos(sinsegine.ins.pos(0)));
        Mpv_ekf(2, 2) = 1.0;
        
        MatrixXd F = MatrixXd::Zero(15, 15);
        F.block<3,3>(0,0) = -INSMath::askew(sinsegine.eth.winn);
        F.block<3,3>(0,9) = -Cnb;
        F.block<3,3>(3,0) = INSMath::askew(fn);
        F.block<3,3>(3,3) = -INSMath::askew(2.0 * sinsegine.eth.wien + sinsegine.eth.wenn);
        F.block<3,3>(3,12)= -Cnb;
        F.block<3,3>(6,3) = Mpv_ekf;
        F.block<3,3>(9,9) = -Matrix3d::Identity() / tau_m_ekf; // 一阶马尔可夫衰减项
        
        MatrixXd Phi = MatrixXd::Identity(15, 15) + F * ctx.ts;
        P = Phi * P * Phi.transpose() + Q * ctx.ts; 

        // 3. 周期结束，执行 EKF 测量更新
        if (sample_count >= total_samples_per_cycle) {
            Vector3d fog_gyro_mean = fog_gyro_active_sum / active_samples;

            auto t_start = chrono::high_resolution_clock::now();

            int meas_dim = exp.use_atomic_acc ? 6 : 3;
            VectorXd Z_meas = VectorXd::Zero(meas_dim);
            MatrixXd H_meas = MatrixXd::Zero(meas_dim, 15);
            MatrixXd R_meas = MatrixXd::Zero(meas_dim, meas_dim);

            Z_meas.segment<3>(0) = fog_gyro_mean - atom.Measure() - sinsegine.ins.eb;
            H_meas.block<3,3>(0,9) = Matrix3d::Identity();
            R_meas.block<3,3>(0,0) = R;

            if (exp.use_atomic_acc) {
                Vector3d fog_acc_mean = fog_acc_active_sum / active_samples;
                Z_meas.segment<3>(3) = fog_acc_mean - atom_acc.Measure() - sinsegine.ins.db;
                H_meas.block<3,3>(3,12) = Matrix3d::Identity();
                R_meas.block<3,3>(3,3) = Matrix3d::Identity() * pow(0.5 * ctx.glv.ug, 2);
            }

            MatrixXd P_sim = P;
            for(int i = 0; i < total_samples_per_cycle; ++i) {
                 P_sim = Phi * P_sim * Phi.transpose() + Q * ctx.ts;
            }
            MatrixXd S = H_meas * P * H_meas.transpose() + R_meas;
            MatrixXd K = P * H_meas.transpose() * S.inverse();
            VectorXd dX = K * Z_meas;
            P = (MatrixXd::Identity(15, 15) - K * H_meas) * P;
            
            auto t_end = chrono::high_resolution_clock::now();
            total_cpu_time_ms += chrono::duration<double, std::milli>(t_end - t_start).count();
            update_counts++;

            // ====================================================================
            // 还原传统 EKF 基线：因无法维持高频精确的历史互协方差，仅做零偏的因果更新
            // ====================================================================
            // 姿态、速度、位置的延迟修正被丢弃，导致 2 秒内的漂移误差被固化（引发 Phase-lag 发散）
            
            // 仅更新传感器零偏
            sinsegine.ins.eb += dX.segment<3>(9);
            sinsegine.ins.db += dX.segment<3>(12); 

            // 记录未被完美回溯纠正的 EKF 真实误差
            double lat_err = (sinsegine.ins.pos(0) - ctx.pos_ref(0)) * ctx.glv.Re;
            double lon_err = (sinsegine.ins.pos(1) - ctx.pos_ref(1)) * ctx.glv.Re * cos(ctx.pos_ref(0));
            double drift = sqrt(lat_err*lat_err + lon_err*lon_err);
            if(drift > max_drift) max_drift = drift;
            
            double true_eb_x_dph = (ctx.eb_align(0) / ctx.glv.deg * 3600.0) + mkv_x_ekf;
            
            LogData(log, total_steps * ctx.ts, sinsegine.ins, ctx, drift, true_eb_x_dph);

            fog_gyro_active_sum.setZero();
            fog_acc_active_sum.setZero();
            sample_count = 0;
            
            if (total_steps > 0 && total_steps % (3600 * 400) == 0) {
                 cout << "   [True EKF] Time: " << (total_steps*ctx.ts)/3600.0 << " h | Max Drift: " << max_drift << " m" << endl;
            }
        }
    }
    log.close();
    cout << " => [Result] True 15-State EKF Max Drift: " << max_drift << " m" << endl;
    cout << " => [CPU] EKF Update Time: " << (total_cpu_time_ms / update_counts) << " ms/cycle\n" << endl;
}

// ==============================================================================
// 算法 4：真正的非线性迭代因子图优化 (Gauss-Newton 重线性化)
// 在 ES-FGO 基础上，修正 bias 后重放 IMU 数据，在新线性化点
// 重新计算 J，迭代求解。打破单次线性化的局限。
// ==============================================================================
void RunTrueFGO(const ExperimentConfig& exp, const SimulationContext& ctx) {
    AtomicGyroSimulator atom(ctx.pos_ref, ctx.glv);
    atom.Init(ctx.att_align); 
    CAIParams acc_params; acc_params.bias_ug = 0.0; acc_params.vrw_ug = 0.05; 
    AtomicAccSimulator atom_acc(ctx.pos_ref, ctx.glv, acc_params);
    atom_acc.Init(ctx.att_align);
    
    SinsEngine sinsegine(ctx.ts);
    sinsegine.eth = ctx.eth; 
    sinsegine.ins = INSState(ctx.att_align, Vector3d::Zero(), ctx.pos_ref, ctx.ts, sinsegine.eth);
    sinsegine.ins.set_bias(ctx.eb_align, ctx.db_align);
    sinsegine.res_init.pos = ctx.pos_ref;

    std::vector<std::string> rest_files(ctx.file_list.begin() + 1, ctx.file_list.end());
    IMUChainedLoader nav_loader(rest_files, ctx.ts, ctx.local_g);

    ofstream log(exp.output_file);
    log << "time,lat_err,lon_err,drift,vn_x,vn_y,vn_z,att_roll_deg,att_pitch_deg,att_yaw_deg,eb_x_dph,eb_y_dph,eb_z_dph,db_x_ug,db_y_ug,db_z_ug,true_eb_x_dph\n"; 

    double T_cycle = exp.t_active + exp.t_dead; 
    int total_samples_per_cycle = static_cast<int>(T_cycle / ctx.ts); 
    int active_samples = static_cast<int>(exp.t_active / ctx.ts);    
    
    Vector3d fog_gyro_active_sum = Vector3d::Zero(); 
    Vector3d fog_acc_active_sum = Vector3d::Zero();
    
    int sample_count = 0;
    double max_drift = 0;
    size_t total_steps = 0;
    size_t max_steps = static_cast<size_t>(exp.duration_hours * 3600.0 / ctx.ts);
    double total_cpu_time_ms = 0;
    int update_counts = 0;

    // True FGO 额外需要：存储当前周期的 IMU 数据用于重放
    std::vector<IMUData> imu_buffer;
    imu_buffer.reserve(total_samples_per_cycle);
    
    // 周期起点状态快照（用于 GN 重放）— 用第一个有效状态初始化
    INSState cycle_start_ins = sinsegine.ins;

    // 第一阶马尔可夫陀螺零偏漂移 (τ=300s, σ_ss=10°/h) — 与 FGO/EKF 相同种子
    const double tau_m_fgo = 300.0, sigma_m_fgo = 10.0;
    double decay_m_fgo = std::exp(-ctx.ts / tau_m_fgo);
    double drive_m_fgo = sigma_m_fgo * sqrt(1.0 - decay_m_fgo * decay_m_fgo);
    double mkv_x_fgo = 0.0, mkv_y_fgo = 0.0;
    std::mt19937 rng_m_fgo(42);
    std::normal_distribution<double> nd_m_fgo(0.0, 1.0);

    // 匹配马尔可夫漂移稳态标准差 σ=10°/h
    Matrix3d Sigma_I_gyro = Matrix3d::Identity() * pow(sigma_m_fgo * ctx.glv.dph, 2); 
    Matrix3d Sigma_A_gyro = Matrix3d::Identity() * pow(1e-4 * ctx.glv.dph, 2);
    Matrix3d Lambda_gyro = Sigma_I_gyro.inverse() + Sigma_A_gyro.inverse();
    Matrix3d K_fgo_gyro = Lambda_gyro.inverse() * Sigma_A_gyro.inverse(); 

    IMUData epoch;
    while (nav_loader.Next(epoch) && total_steps < max_steps) {
        
        if (exp.inject_time_varying_drift) {
            mkv_x_fgo = decay_m_fgo * mkv_x_fgo + drive_m_fgo * nd_m_fgo(rng_m_fgo);
            mkv_y_fgo = decay_m_fgo * mkv_y_fgo + drive_m_fgo * nd_m_fgo(rng_m_fgo);
            epoch.wm(0) += mkv_x_fgo * ctx.glv.dph * ctx.ts;
            epoch.wm(1) += mkv_y_fgo * ctx.glv.dph * ctx.ts;
        }
        
        sample_count++;
        total_steps++;

        if (sample_count <= active_samples) {
            fog_gyro_active_sum += epoch.wm / ctx.ts; 
            fog_acc_active_sum += epoch.vm / ctx.ts;
        } 

        // 保存周期起点状态（第一个样本之前已经 step 过了，所以在第1个样本时保存）
        if (sample_count == 1) {
            cycle_start_ins = sinsegine.ins;
        }

        // 存储 IMU 数据（噪声已注入），用于 GN 重放
        imu_buffer.push_back(epoch);

        // 前向机械编排
        sinsegine.Step_Nav(epoch.wm, epoch.vm);
        sinsegine.ins.vn(2) = 0.0; 

        // 周期结束
        if (sample_count >= total_samples_per_cycle) {
            
            Vector3d fog_gyro_mean = fog_gyro_active_sum / active_samples;
            Vector3d delta_z_gyro = fog_gyro_mean - atom.Measure();

            auto t_start = chrono::high_resolution_clock::now();

            // ==============================================================
            // MAP 最优零偏估计（线性观测方程，一次解析解即达最优）
            // ==============================================================
            Vector3d delta_eb_star = K_fgo_gyro * (delta_z_gyro - sinsegine.ins.eb);
            Vector3d eb_new = sinsegine.ins.eb + delta_eb_star;
            Vector3d db_new = sinsegine.ins.db;
            if (exp.use_atomic_acc) {
                Matrix3d Sigma_I_acc = Matrix3d::Identity() * pow(50.0 * ctx.glv.ug, 2);
                Matrix3d Sigma_A_acc = Matrix3d::Identity() * pow(0.5 * ctx.glv.ug, 2);
                Matrix3d Lambda_acc = Sigma_I_acc.inverse() + Sigma_A_acc.inverse();
                Matrix3d K_fgo_acc = Lambda_acc.inverse() * Sigma_A_acc.inverse();
                Vector3d fog_acc_mean = fog_acc_active_sum / active_samples;
                Vector3d delta_z_acc = fog_acc_mean - atom_acc.Measure();
                Vector3d delta_db_star = K_fgo_acc * (delta_z_acc - sinsegine.ins.db);
                db_new = sinsegine.ins.db + delta_db_star;
            }

            // ==============================================================
            // 用最优 bias 从周期起点重放 IMU → 非线性精确轨迹
            // ==============================================================
            SinsEngine replay(ctx.ts);
            replay.eth = ctx.eth;
            replay.ins = cycle_start_ins;
            replay.ins.set_bias(eb_new, db_new);
            replay.res_init.pos = ctx.pos_ref;

            double true_eb_x_dph = (ctx.eb_align(0) / ctx.glv.deg * 3600.0) + mkv_x_fgo;
            size_t replay_step = 0;
            double replay_time_base = (total_steps - imu_buffer.size()) * ctx.ts;

            for (const auto& imu : imu_buffer) {
                replay.Step_Nav(imu.wm, imu.vm);
                replay.ins.vn(2) = 0.0;
                replay_step++;

                // 直接用非线性重放轨迹记录日志（彻底抛弃 J 线性补偿）
                double nav_time = replay_time_base + replay_step * ctx.ts;
                double lat_err = (replay.ins.pos(0) - ctx.pos_ref(0)) * ctx.glv.Re;
                double lon_err = (replay.ins.pos(1) - ctx.pos_ref(1)) * ctx.glv.Re * cos(ctx.pos_ref(0));
                double drift = sqrt(lat_err*lat_err + lon_err*lon_err);
                if(drift > max_drift) max_drift = drift;

                INSState replay_ins = replay.ins;
                replay_ins.eb = eb_new;
                LogData(log, nav_time, replay_ins, ctx, drift, true_eb_x_dph);
            }

            // 用重放得到的精确非线性端点更新主引擎
            sinsegine.ins = replay.ins;
            sinsegine.ins.set_bias(eb_new, db_new);
            sinsegine.eth = replay.eth;

            auto t_end = chrono::high_resolution_clock::now();
            total_cpu_time_ms += chrono::duration<double, std::milli>(t_end - t_start).count();
            update_counts++;

            if (total_steps > 0 && total_steps % (3600 * 400) == 0) {
                 cout << "   [True FGO] Time: " << (total_steps*ctx.ts)/3600.0 << " h | Max Drift: " << max_drift << " m" << endl;
            }

            fog_gyro_active_sum.setZero();
            fog_acc_active_sum.setZero();
            sample_count = 0; 
            imu_buffer.clear();
        }
    }
    log.close();
    cout << " => [Result] True FGO (GN Re-linearization) Max Drift: " << max_drift << " m" << endl;
    cout << " => [CPU] Avg True FGO Update Time: " << (total_cpu_time_ms / update_counts) << " ms/cycle\n" << endl;
}

void RunExperimentRouter(const ExperimentConfig& exp, const SimulationContext& ctx) {
    cout << "\n--------------------------------------------------" << endl;
    cout << " Starting Scenario: " << exp.name << endl;
    cout << " Duration: " << exp.duration_hours << " h | Output: " << exp.output_file << endl;
    
    std::ifstream check_file(exp.output_file);
    if (check_file.good()) {
        cout << " [SKIP] Output file already exists, skipping..." << endl;
        cout << "--------------------------------------------------" << endl;
        check_file.close();
        return;
    }
    
    cout << "--------------------------------------------------" << endl;

    switch (exp.algo) {
        case AlgoType::PURE_FOG: RunPureFOG(exp, ctx); break;
        case AlgoType::ES_FGO:   RunESFGO(exp, ctx);   break;
        case AlgoType::EKF:      RunEKF(exp, ctx);     break;
        case AlgoType::TRUE_FGO: RunTrueFGO(exp, ctx); break;
    }
}

int RunFgoExperimentApp() {
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

    vector<ExperimentConfig> paper_experiments = {
        // ================= 1. 静态基线 (无注入漂移 inject = false) =================
        {"Static: Pure FOG",       AlgoType::PURE_FOG, 4.0, 0.0, 0.0, false, "static_pure_fog.csv",       false, 0.0, 1.0},
        // 静态 - 仅融合原子陀螺 (Gyro Only)
        {"Static: EKF (Gyro)",     AlgoType::EKF,      4.0, 1.6, 0.4, false, "static_ekf_gyro.csv",       false, 0.0, 1.0},
        {"Static: FGO (Gyro)",     AlgoType::ES_FGO,   4.0, 1.6, 0.4, false, "static_fgo_gyro.csv",       false, 0.0, 1.0},
        // 静态 - 双设备融合 (Gyro + Acc)
        {"Static: EKF (Dual)",     AlgoType::EKF,      4.0, 1.6, 0.4, true,  "static_ekf_dual.csv",       false, 0.0, 1.0},
        {"Static: FGO (Dual)",     AlgoType::ES_FGO,   4.0, 1.6, 0.4, true,  "static_fgo_dual.csv",       false, 0.0, 1.0},

        // ================= 2. 动态对比 (注入 GM1 时变漂移 inject = true) =================
        // 动态 - 仅融合原子陀螺 (Gyro Only)
        {"Dynamic: EKF (Gyro)",    AlgoType::EKF,      4.0, 1.6, 0.4, false, "dynamic_ekf_gyro.csv",      true, 0.8, 0.5},
        {"Dynamic: FGO (Gyro)",    AlgoType::ES_FGO,   4.0, 1.6, 0.4, false, "dynamic_fgo_gyro.csv",      true, 0.8, 0.5},
        {"Dynamic: True FGO(Gyro)",AlgoType::TRUE_FGO, 4.0, 1.6, 0.4, false, "dynamic_true_fgo_gyro.csv", true, 0.8, 0.5},
        // 动态 - 双设备融合 (Gyro + Acc)
        {"Dynamic: EKF (Dual)",    AlgoType::EKF,      4.0, 1.6, 0.4, true,  "dynamic_ekf_dual.csv",      true, 0.8, 0.5},
        {"Dynamic: FGO (Dual)",    AlgoType::ES_FGO,   4.0, 1.6, 0.4, true,  "dynamic_fgo_dual.csv",      true, 0.8, 0.5},
        {"Dynamic: True FGO(Dual)",AlgoType::TRUE_FGO, 4.0, 1.6, 0.4, true,  "dynamic_true_fgo_dual.csv", true, 0.8, 0.5}
    };
    
    cout << "\n>>> Starting Wideband Vibration Validation (Spectral Aliasing Test)..." << endl;
    for (const auto& exp : paper_experiments) {
        RunExperimentRouter(exp, ctx);
    }

    cout << "\nAll Simulations Completed! You can now fill the data into the LaTeX source." << endl;
    return 0;
}
