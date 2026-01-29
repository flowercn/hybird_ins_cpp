/**
 * @file fgo_sim_test.cpp
 * @brief 仿真数据验证因子图优化理论极限
 * 
 * 用两组仿真数据对比：
 *   1. 纯ARW陀螺 - 验证理论公式
 *   2. ARW + 常值零偏 - 验证零偏可观测性
 */

#include <iostream>
#include <iomanip>
#include <random>
#include <vector>
#include <chrono>

// GTSAM
#include <gtsam/navigation/CombinedImuFactor.h>
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
using gtsam::symbol_shorthand::X;
using gtsam::symbol_shorthand::V;
using gtsam::symbol_shorthand::B;

constexpr double kDeg2Rad = M_PI / 180.0;
constexpr double kRad2Deg = 180.0 / M_PI;

// ============ 生成仿真 IMU 数据 ============
struct SimConfig {
    double ts = 1.0 / 400.0;
    double duration = 3600.0;  // 1小时
    double arw_deg_sqrt_h = 0.003;  // deg/sqrt(h)
    double bias_deg_h[3] = {0, 0, 0};  // deg/h
    double g = 9.8;
    double lat = 32.0286 * kDeg2Rad;
};

vector<IMUData> GenerateSimIMU(const SimConfig& cfg, mt19937& rng) {
    size_t n = static_cast<size_t>(cfg.duration / cfg.ts);
    vector<IMUData> data(n);
    
    // ARW -> 离散噪声标准差
    // ARW (deg/√h) -> σ_ω (rad/s) = ARW * π/180 / √(3600 * fs)
    double fs = 1.0 / cfg.ts;
    double sigma_w = cfg.arw_deg_sqrt_h * kDeg2Rad / sqrt(3600.0 * fs);  // rad/s
    
    // 地球自转角速度（导航系）
    double wie = 7.2921151467e-5;  // rad/s
    Vector3d wie_n(0, wie * cos(cfg.lat), wie * sin(cfg.lat));
    
    // 重力在导航系（假设北东地）
    Vector3d gn(0, 0, cfg.g);
    
    // 零偏
    Vector3d eb(cfg.bias_deg_h[0] * kDeg2Rad / 3600.0,
                cfg.bias_deg_h[1] * kDeg2Rad / 3600.0,
                cfg.bias_deg_h[2] * kDeg2Rad / 3600.0);
    
    normal_distribution<double> noise_w(0.0, sigma_w);
    normal_distribution<double> noise_a(0.0, 50e-6 * cfg.g);  // 50 ug 白噪声
    
    for (size_t i = 0; i < n; ++i) {
        // 陀螺输出 = 地球自转 + 零偏 + 白噪声
        // 注意：体坐标系和导航系对齐时 wie_b = wie_n
        Vector3d w = wie_n + eb;
        w(0) += noise_w(rng);
        w(1) += noise_w(rng);
        w(2) += noise_w(rng);
        
        // 角增量
        data[i].wm = w * cfg.ts;
        
        // 加速度计：静止时测量比力 = -g (在导航系向下为正时)
        // 但北东地坐标系下，z轴向下，所以比力向上 = (0, 0, g)
        Vector3d f = gn;  // 比力
        f(0) += noise_a(rng);
        f(1) += noise_a(rng);
        f(2) += noise_a(rng);
        data[i].vm = f * cfg.ts;
    }
    
    return data;
}

// ============ 传统惯导解算 ============
double RunTraditionalINS(const vector<IMUData>& imu_data, const SimConfig& cfg,
                         const Vector3d& att0, const Vector3d& eb_est) {
    GLV glv;
    Vector3d pos0(cfg.lat, 118.8533 * glv.deg, 17.0);
    
    SinsEngine ins(cfg.ts);
    ins.Sins_Init(pos0, Vector3d::Zero(), att0, eb_est, Vector3d::Zero());
    
    double max_drift = 0;
    for (const auto& imu : imu_data) {
        ins.Step_Nav(imu.wm, imu.vm);
        
        double dN = (ins.ins.pos(0) - pos0(0)) * glv.Re;
        double dE = (ins.ins.pos(1) - pos0(1)) * cos(pos0(0)) * glv.Re;
        double drift = sqrt(dN*dN + dE*dE);
        if (drift > max_drift) max_drift = drift;
    }
    
    return max_drift / 1852.0;  // nm
}

// ============ FGO 优化 ============
double RunFGO(const vector<IMUData>& imu_data, const SimConfig& cfg,
              const Vector3d& att0, const Vector3d& eb_est) {
    
    GLV glv;
    Vector3d pos0(cfg.lat, 118.8533 * glv.deg, 17.0);
    
    // 预积分参数
    double arw_rad = cfg.arw_deg_sqrt_h * kDeg2Rad / 60.0;  // rad/s/sqrt(Hz)
    auto preint_params = PreintegrationCombinedParams::MakeSharedU(cfg.g);
    preint_params->setGyroscopeCovariance(pow(arw_rad, 2) * Eigen::Matrix3d::Identity());
    preint_params->setAccelerometerCovariance(1e-6 * Eigen::Matrix3d::Identity());
    preint_params->setIntegrationCovariance(1e-8 * Eigen::Matrix3d::Identity());
    preint_params->setBiasAccCovariance(1e-10 * Eigen::Matrix3d::Identity());
    preint_params->setBiasOmegaCovariance(1e-10 * Eigen::Matrix3d::Identity());
    preint_params->setBiasAccOmegaInit(1e-10 * Eigen::Matrix<double,6,6>::Identity());
    
    imuBias::ConstantBias prior_bias(Vector3d::Zero(), eb_est);
    
    Rot3 R0 = Rot3::RzRyRx(att0(0), att0(1), att0(2));
    Pose3 pose0(R0, Point3(pos0(0), pos0(1), pos0(2)));
    Vector3d vel0 = Vector3d::Zero();
    Vector3d zero_vel = Vector3d::Zero();
    
    NonlinearFactorGraph graph;
    Values initial_values;
    
    // 先验
    auto noise_pose = noiseModel::Diagonal::Sigmas(
        (Vector6() << 0.001, 0.001, 0.005, 0.001, 0.001, 0.001).finished());
    auto noise_vel = noiseModel::Isotropic::Sigma(3, 0.01);
    auto noise_bias = noiseModel::Diagonal::Sigmas(
        (Vector6() << 1e-5, 1e-5, 1e-5, 1e-7, 1e-7, 1e-7).finished());
    
    graph.addPrior(X(0), pose0, noise_pose);
    graph.addPrior(V(0), vel0, noise_vel);
    graph.addPrior(B(0), prior_bias, noise_bias);
    
    initial_values.insert(X(0), pose0);
    initial_values.insert(V(0), vel0);
    initial_values.insert(B(0), prior_bias);
    
    auto noise_zero_vel = noiseModel::Isotropic::Sigma(3, 0.001);
    
    // 每秒一个节点
    size_t node_samples = static_cast<size_t>(1.0 / cfg.ts);
    auto preint = std::make_shared<PreintegratedCombinedMeasurements>(preint_params, prior_bias);
    
    size_t node_idx = 0;
    for (size_t i = 0; i < imu_data.size(); ++i) {
        Vector3d omega = imu_data[i].wm / cfg.ts;
        Vector3d accel = imu_data[i].vm / cfg.ts;
        preint->integrateMeasurement(accel, omega, cfg.ts);
        
        if ((i + 1) % node_samples == 0) {
            size_t next_node = node_idx + 1;
            
            graph.emplace_shared<CombinedImuFactor>(
                X(node_idx), V(node_idx), X(next_node), V(next_node), B(node_idx), B(next_node),
                *preint);
            
            graph.addPrior(V(next_node), zero_vel, noise_zero_vel);
            
            // 零位移约束
            auto noise_pos_only = noiseModel::Diagonal::Sigmas(
                (Vector6() << 1e3, 1e3, 1e3, 0.01, 0.01, 0.01).finished());
            graph.addPrior(X(next_node), pose0, noise_pos_only);
            
            NavState prev_state(pose0, vel0);
            NavState pred_state = preint->predict(prev_state, prior_bias);
            
            initial_values.insert(X(next_node), pred_state.pose());
            initial_values.insert(V(next_node), pred_state.velocity());
            initial_values.insert(B(next_node), prior_bias);
            
            preint->resetIntegrationAndSetBias(prior_bias);
            node_idx = next_node;
        }
    }
    
    // 优化
    LevenbergMarquardtParams opt_params;
    opt_params.setVerbosity("SILENT");
    opt_params.setMaxIterations(30);
    
    LevenbergMarquardtOptimizer optimizer(graph, initial_values, opt_params);
    Values result = optimizer.optimize();
    
    // 计算漂移
    Point3 pos_init(pos0(0), pos0(1), pos0(2));
    double max_drift = 0;
    
    for (size_t k = 0; k <= node_idx; ++k) {
        Pose3 pose_k = result.at<Pose3>(X(k));
        double dN = (pose_k.translation().x() - pos_init.x()) * glv.Re;
        double dE = (pose_k.translation().y() - pos_init.y()) * cos(cfg.lat) * glv.Re;
        double drift = sqrt(dN*dN + dE*dE);
        if (drift > max_drift) max_drift = drift;
    }
    
    return max_drift / 1852.0;  // nm
}

int main() {
    cout << "============================================================\n";
    cout << "  仿真验证：因子图优化理论极限\n";
    cout << "============================================================\n";
    
    mt19937 rng(42);  // 固定种子
    SimConfig cfg;
    
    // ========== Case 1: 纯 ARW，无零偏 ==========
    cout << "\n[Case 1] 纯 ARW = " << cfg.arw_deg_sqrt_h << " deg/√h, 无零偏\n";
    
    auto imu1 = GenerateSimIMU(cfg, rng);
    Vector3d att0 = Vector3d::Zero();  // 完美对准
    Vector3d eb_est = Vector3d::Zero();  // 完美零偏估计
    
    double trad1 = RunTraditionalINS(imu1, cfg, att0, eb_est);
    double fgo1 = RunFGO(imu1, cfg, att0, eb_est);
    
    // 理论值
    double arw_rad = cfg.arw_deg_sqrt_h * kDeg2Rad;
    double theory1 = (1.0/6.0) * cfg.g * arw_rad * pow(cfg.duration, 2) / 1852.0;
    
    cout << "    理论极限: " << fixed << setprecision(4) << theory1 << " nm\n";
    cout << "    传统INS:  " << trad1 << " nm\n";
    cout << "    FGO优化:  " << fgo1 << " nm\n";
    
    // ========== Case 2: ARW + 未知零偏 ==========
    cout << "\n[Case 2] ARW + 未知零偏 0.01 deg/h\n";
    
    SimConfig cfg2 = cfg;
    cfg2.bias_deg_h[0] = 0.01;
    cfg2.bias_deg_h[1] = 0.01;
    cfg2.bias_deg_h[2] = 0.01;
    
    auto imu2 = GenerateSimIMU(cfg2, rng);
    
    // 传统方法：不知道零偏
    double trad2_no_eb = RunTraditionalINS(imu2, cfg2, att0, Vector3d::Zero());
    // 传统方法：知道零偏
    Vector3d true_eb(0.01 * kDeg2Rad / 3600, 0.01 * kDeg2Rad / 3600, 0.01 * kDeg2Rad / 3600);
    double trad2_with_eb = RunTraditionalINS(imu2, cfg2, att0, true_eb);
    // FGO：不给零偏
    double fgo2 = RunFGO(imu2, cfg2, att0, Vector3d::Zero());
    
    cout << "    传统INS (无零偏补偿):  " << fixed << setprecision(4) << trad2_no_eb << " nm\n";
    cout << "    传统INS (有零偏补偿):  " << trad2_with_eb << " nm\n";
    cout << "    FGO优化 (估计零偏):    " << fgo2 << " nm\n";
    
    // ========== Case 3: 初始姿态误差 ==========
    cout << "\n[Case 3] 初始姿态误差 0.01 deg (无零偏)\n";
    
    auto imu3 = GenerateSimIMU(cfg, rng);  // 纯ARW数据
    Vector3d att_err(0.01 * kDeg2Rad, 0.01 * kDeg2Rad, 0);  // 水平姿态误差
    
    double trad3 = RunTraditionalINS(imu3, cfg, att_err, Vector3d::Zero());
    double fgo3 = RunFGO(imu3, cfg, att_err, Vector3d::Zero());
    
    // 理论：常值姿态误差
    double theta = 0.01 * kDeg2Rad * sqrt(2);  // 水平合成
    double theory3 = 0.5 * cfg.g * theta * pow(cfg.duration, 2) / 1852.0;
    
    cout << "    理论 (常值姿态误差): " << fixed << setprecision(4) << theory3 << " nm\n";
    cout << "    传统INS:             " << trad3 << " nm\n";
    cout << "    FGO优化:             " << fgo3 << " nm\n";
    
    // ========== 统计 Monte Carlo ==========
    cout << "\n[Case 4] Monte Carlo 统计 (20次, 纯ARW)\n";
    
    vector<double> trad_results, fgo_results;
    for (int trial = 0; trial < 20; ++trial) {
        auto imu_mc = GenerateSimIMU(cfg, rng);
        trad_results.push_back(RunTraditionalINS(imu_mc, cfg, Vector3d::Zero(), Vector3d::Zero()));
        fgo_results.push_back(RunFGO(imu_mc, cfg, Vector3d::Zero(), Vector3d::Zero()));
    }
    
    auto mean = [](const vector<double>& v) {
        return accumulate(v.begin(), v.end(), 0.0) / v.size();
    };
    auto stddev = [&mean](const vector<double>& v) {
        double m = mean(v);
        double sq_sum = 0;
        for (auto x : v) sq_sum += (x - m) * (x - m);
        return sqrt(sq_sum / v.size());
    };
    
    cout << "    传统INS: " << fixed << setprecision(4) 
         << mean(trad_results) << " ± " << stddev(trad_results) << " nm\n";
    cout << "    FGO优化: " 
         << mean(fgo_results) << " ± " << stddev(fgo_results) << " nm\n";
    cout << "    理论:    " << theory1 << " nm\n";
    
    cout << "\n============================================================\n";
    
    return 0;
}
