/**
 * @file fgo_test.cpp
 * @brief 因子图优化(FGO)惯导测试
 * 
 * 核心思想：
 *   - 使用 GTSAM 的 IMU 预积分因子
 *   - 静止状态：零速约束 + 零位移约束
 *   - 批量优化（非增量）
 * 
 * 因子图结构：
 *   X0 ---- preint ---- X1 ---- preint ---- X2 ...
 *   |                   |                   |
 *   prior               zero_vel            zero_vel
 *                       zero_pos            zero_pos
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>

// GTSAM
#include <gtsam/navigation/CombinedImuFactor.h>
#include <gtsam/navigation/ImuFactor.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/slam/PriorFactor.h>

// Local
#include "sins_engine.h"
#include "support.h"

using namespace std;
using namespace Eigen;
using namespace gtsam;
using gtsam::symbol_shorthand::X;  // Pose3 (position, rotation)
using gtsam::symbol_shorthand::V;  // Vector3 (velocity)
using gtsam::symbol_shorthand::B;  // imuBias::ConstantBias

// ============ 配置 ============
constexpr double kDeg2Rad = M_PI / 180.0;
constexpr double kRad2Deg = 180.0 / M_PI;

struct FGOConfig {
    double ts = 1.0 / 400.0;          // IMU 采样周期
    double node_interval = 1.0;       // 因子图节点间隔 (秒)
    double align_time = 3600.0;       // 对准时间 (秒)
    double nav_time = 3600.0;         // 导航时间 (秒)
    
    // IMU 噪声参数 (连续时间)
    double gyro_noise_density = 0.003 * kDeg2Rad / 60.0;    // rad/s/sqrt(Hz) ~ 0.003 deg/sqrt(h)
    double accel_noise_density = 50e-6 * 9.8;               // m/s^2/sqrt(Hz) ~ 50 ug/sqrt(Hz)
    double gyro_bias_rw = 0.0003 * kDeg2Rad / 3600.0;       // rad/s^2/sqrt(Hz)
    double accel_bias_rw = 1e-6 * 9.8;                      // m/s^3/sqrt(Hz)
    
    // 约束噪声
    double zero_vel_sigma = 0.001;    // m/s
    double zero_pos_sigma = 0.01;     // m
};

// ============ 计算漂移 ============
double ComputeDrift(const Point3& pos, const Point3& pos0, const GLV& glv) {
    double lat0 = pos0.x();
    double dN = (pos.x() - pos0.x()) * glv.Re;
    double dE = (pos.y() - pos0.y()) * cos(lat0) * glv.Re;
    return sqrt(dN*dN + dE*dE);
}

int main() {
    auto t_start = chrono::high_resolution_clock::now();
    
    cout << "=============================================================\n";
    cout << "  Factor Graph Optimization (FGO) INS Test\n";
    cout << "=============================================================\n\n";
    
    FGOConfig cfg;
    GLV glv;
    
    // ========== 1. 加载数据 ==========
    cout << "[1] Loading data..." << endl;
    
    Vector3d pos0(32.0286 * glv.deg, 118.8533 * glv.deg, 17.0);
    Earth eth(glv);
    eth.update(pos0, Vector3d::Zero());
    double local_g = eth.gn.norm();
    
    auto imu_data = LoadIMUData("../fog3h.csv", cfg.ts, local_g);
    cout << "    Total: " << imu_data.size() << " samples ("
         << imu_data.size() * cfg.ts / 3600 << " h)\n";
    
    // 应用刻度因数修正
    double K_scale = 1.000372;
    for (auto& imu : imu_data) {
        imu.wm *= K_scale;
    }
    cout << "    Scale factor K = " << K_scale << " applied\n";
    
    // ========== 2. 传统方法对准（获取初始姿态）==========
    cout << "\n[2] Coarse alignment for initial attitude..." << endl;
    
    SinsEngine ins(cfg.ts);
    ins.Sins_Init(pos0, Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero());
    
    HybridAlignConfig align_cfg;
    align_cfg.t_coarse = 60;
    align_cfg.t_fine = cfg.align_time;
    align_cfg.eb_sigma_allan = 0.003;
    align_cfg.verbose = false;
    
    auto align_result = ins.Run_HybridAlign(imu_data, align_cfg);
    if (!align_result.valid) {
        cerr << "Alignment failed!" << endl;
        return -1;
    }
    
    Vector3d att0 = align_result.att;
    cout << "    Initial attitude (deg): " 
         << att0(0)/glv.deg << ", " << att0(1)/glv.deg << ", " << att0(2)/glv.deg << endl;
    
    // ========== 3. 构建因子图 ==========
    cout << "\n[3] Building factor graph..." << endl;
    
    // 预积分参数
    auto preint_params = PreintegrationCombinedParams::MakeSharedU(local_g);
    preint_params->setGyroscopeCovariance(pow(cfg.gyro_noise_density, 2) * Eigen::Matrix3d::Identity());
    preint_params->setAccelerometerCovariance(pow(cfg.accel_noise_density, 2) * Eigen::Matrix3d::Identity());
    preint_params->setIntegrationCovariance(1e-8 * Eigen::Matrix3d::Identity());
    preint_params->setBiasAccCovariance(pow(cfg.accel_bias_rw, 2) * Eigen::Matrix3d::Identity());
    preint_params->setBiasOmegaCovariance(pow(cfg.gyro_bias_rw, 2) * Eigen::Matrix3d::Identity());
    preint_params->setBiasAccOmegaInit(1e-5 * Eigen::Matrix<double,6,6>::Identity());
    
    // 初始偏差
    imuBias::ConstantBias prior_bias(Vector3d::Zero(), align_result.eb);
    
    // 初始姿态（注意：GTSAM 的 Rot3::RzRyRx 顺序是 roll, pitch, yaw）
    Rot3 R0 = Rot3::RzRyRx(att0(0), att0(1), att0(2));  // pitch, roll, yaw
    Pose3 pose0(R0, Point3(pos0(0), pos0(1), pos0(2)));
    Vector3d vel0 = Vector3d::Zero();
    
    NonlinearFactorGraph graph;
    Values initial_values;
    
    // 先验因子
    auto noise_pose = noiseModel::Diagonal::Sigmas(
        (Vector6() << 0.001, 0.001, 0.005, 0.001, 0.001, 0.001).finished());  // rad, rad, rad, m, m, m
    auto noise_vel = noiseModel::Isotropic::Sigma(3, 0.01);
    auto noise_bias = noiseModel::Diagonal::Sigmas(
        (Vector6() << 1e-5, 1e-5, 1e-5, 1e-7, 1e-7, 1e-7).finished());  // accel, gyro
    
    graph.addPrior(X(0), pose0, noise_pose);
    graph.addPrior(V(0), vel0, noise_vel);
    graph.addPrior(B(0), prior_bias, noise_bias);
    
    initial_values.insert(X(0), pose0);
    initial_values.insert(V(0), vel0);
    initial_values.insert(B(0), prior_bias);
    
    // 零速零位移约束噪声
    auto noise_zero_vel = noiseModel::Isotropic::Sigma(3, cfg.zero_vel_sigma);
    auto noise_zero_pos = noiseModel::Isotropic::Sigma(3, cfg.zero_pos_sigma);
    
    // 导航起始索引
    size_t nav_start_idx = static_cast<size_t>(cfg.align_time / cfg.ts);
    size_t nav_end_idx = nav_start_idx + static_cast<size_t>(cfg.nav_time / cfg.ts);
    if (nav_end_idx > imu_data.size()) nav_end_idx = imu_data.size();
    
    // 节点间隔（采样数）
    size_t node_samples = static_cast<size_t>(cfg.node_interval / cfg.ts);
    
    // 预积分并添加因子
    auto preint = std::make_shared<PreintegratedCombinedMeasurements>(preint_params, prior_bias);
    
    size_t node_idx = 0;
    size_t num_nodes = (nav_end_idx - nav_start_idx) / node_samples;
    
    // 零速零位移的先验值
    Vector3d zero_vel = Vector3d::Zero();
    
    cout << "    Nodes: " << num_nodes << " (interval: " << cfg.node_interval << " s)\n";
    
    for (size_t i = nav_start_idx; i < nav_end_idx; ++i) {
        // 获取 IMU 数据（角增量/速度增量 -> 角速度/加速度）
        Vector3d omega = imu_data[i].wm / cfg.ts;
        Vector3d accel = imu_data[i].vm / cfg.ts;
        
        preint->integrateMeasurement(accel, omega, cfg.ts);
        
        // 到达节点时刻
        if ((i - nav_start_idx + 1) % node_samples == 0) {
            size_t next_node = node_idx + 1;
            
            // 添加 IMU 预积分因子
            graph.emplace_shared<CombinedImuFactor>(
                X(node_idx), V(node_idx), X(next_node), V(next_node), B(node_idx), B(next_node),
                *preint);
            
            // 添加零速约束
            graph.addPrior(V(next_node), zero_vel, noise_zero_vel);
            
            // 添加零位移约束（位置应该不变）
            // 通过对 Pose 施加位置先验实现
            // 使用 Pose 先验但位置噪声小、姿态噪声大
            auto noise_pos_only = noiseModel::Diagonal::Sigmas(
                (Vector6() << 1e3, 1e3, 1e3, cfg.zero_pos_sigma, cfg.zero_pos_sigma, cfg.zero_pos_sigma).finished());
            graph.addPrior(X(next_node), pose0, noise_pos_only);
            
            // 预测下一状态作为初始值
            NavState prev_state(pose0, vel0);  // 静止假设
            NavState pred_state = preint->predict(prev_state, prior_bias);
            
            initial_values.insert(X(next_node), pred_state.pose());
            initial_values.insert(V(next_node), pred_state.velocity());
            initial_values.insert(B(next_node), prior_bias);
            
            // 重置预积分
            preint->resetIntegrationAndSetBias(prior_bias);
            node_idx = next_node;
            
            if (node_idx % 100 == 0) {
                cout << "    Built " << node_idx << " / " << num_nodes << " nodes\r" << flush;
            }
        }
    }
    cout << "    Built " << node_idx << " nodes total              \n";
    
    // ========== 4. 优化 ==========
    cout << "\n[4] Optimizing..." << endl;
    
    LevenbergMarquardtParams opt_params;
    opt_params.setVerbosityLM("SUMMARY");
    opt_params.setMaxIterations(50);
    
    auto t_opt_start = chrono::high_resolution_clock::now();
    
    LevenbergMarquardtOptimizer optimizer(graph, initial_values, opt_params);
    Values result = optimizer.optimize();
    
    auto t_opt_end = chrono::high_resolution_clock::now();
    double opt_time = chrono::duration<double>(t_opt_end - t_opt_start).count();
    
    cout << "    Optimization time: " << fixed << setprecision(2) << opt_time << " s\n";
    cout << "    Final error: " << graph.error(result) << endl;
    
    // ========== 5. 分析结果 ==========
    cout << "\n[5] Results:" << endl;
    
    // 计算位置漂移
    Point3 pos_init(pos0(0), pos0(1), pos0(2));
    double max_drift = 0;
    
    cout << "    Time(min)  Drift(m)  Vel(m/s)\n";
    for (size_t k = 0; k <= node_idx; k += max(1UL, node_idx / 10)) {
        Pose3 pose_k = result.at<Pose3>(X(k));
        Vector3d vel_k = result.at<Vector3d>(V(k));
        
        double drift = ComputeDrift(pose_k.translation(), pos_init, glv);
        if (drift > max_drift) max_drift = drift;
        
        double t_min = k * cfg.node_interval / 60.0;
        cout << "    " << fixed << setw(8) << setprecision(1) << t_min
             << setw(10) << setprecision(3) << drift
             << setw(10) << setprecision(4) << vel_k.norm() << endl;
    }
    
    // 最终结果
    Pose3 final_pose = result.at<Pose3>(X(node_idx));
    imuBias::ConstantBias final_bias = result.at<imuBias::ConstantBias>(B(node_idx));
    
    double final_drift = ComputeDrift(final_pose.translation(), pos_init, glv);
    double drift_rate = max_drift / (cfg.nav_time / 3600.0) / 1852.0;
    
    cout << "\n    Final drift: " << fixed << setprecision(3) << final_drift << " m\n";
    cout << "    Max drift:   " << max_drift << " m\n";
    cout << "    Drift rate:  " << setprecision(4) << drift_rate << " nm/h\n";
    
    // 估计的零偏
    Vector3d gyro_bias = final_bias.gyroscope();
    Vector3d accel_bias = final_bias.accelerometer();
    cout << "\n    Estimated gyro bias (deg/h): "
         << gyro_bias(0) * kRad2Deg * 3600 << ", "
         << gyro_bias(1) * kRad2Deg * 3600 << ", "
         << gyro_bias(2) * kRad2Deg * 3600 << endl;
    cout << "    Estimated accel bias (ug): "
         << accel_bias(0) / 9.8 * 1e6 << ", "
         << accel_bias(1) / 9.8 * 1e6 << ", "
         << accel_bias(2) / 9.8 * 1e6 << endl;
    
    // 总时间
    auto t_end = chrono::high_resolution_clock::now();
    double total_time = chrono::duration<double>(t_end - t_start).count();
    cout << "\n    Total time: " << fixed << setprecision(1) << total_time << " s\n";
    
    // ========== 6. 对比传统方法 ==========
    cout << "\n[6] Comparison with traditional INS..." << endl;
    
    ins.eth.update(pos0, Vector3d::Zero());
    ins.ins = INSState(align_result.att, Vector3d::Zero(), pos0, cfg.ts, ins.eth);
    ins.ins.set_bias(align_result.eb, align_result.db);
    
    double trad_max_drift = 0;
    for (size_t i = nav_start_idx; i < nav_end_idx; ++i) {
        ins.Step_Nav(imu_data[i].wm, imu_data[i].vm);
        
        double dN = (ins.ins.pos(0) - pos0(0)) * glv.Re;
        double dE = (ins.ins.pos(1) - pos0(1)) * cos(pos0(0)) * glv.Re;
        double drift = sqrt(dN*dN + dE*dE);
        if (drift > trad_max_drift) trad_max_drift = drift;
    }
    double trad_drift_rate = trad_max_drift / (cfg.nav_time / 3600.0) / 1852.0;
    
    cout << "    Traditional INS:  " << fixed << setprecision(4) << trad_drift_rate << " nm/h\n";
    cout << "    FGO Optimized:    " << drift_rate << " nm/h\n";
    
    if (drift_rate < trad_drift_rate) {
        double improve = (trad_drift_rate - drift_rate) / trad_drift_rate * 100;
        cout << "    Improvement:      " << setprecision(1) << improve << "%\n";
    } else {
        cout << "    (No improvement)\n";
    }
    
    cout << "\n=============================================================\n";
    
    return 0;
}
