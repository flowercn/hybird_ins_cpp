#include "fgo_core.h"
#include <iostream>
#include <boost/make_shared.hpp>
#include <gtsam/linear/linearExceptions.h> 

using namespace std;
using namespace gtsam;
using namespace gtsam::symbol_shorthand;

FGOEngine::FGOEngine(const FGOSettings& settings) : settings_(settings) {
    Eigen::Vector3d n_gravity(0, 0, -settings_.gravity); 
    p_params_ = boost::make_shared<PreintegratedCombinedMeasurements::Params>(n_gravity);
    
    Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
    p_params_->gyroscopeCovariance = I3 * pow(settings_.gyro_noise, 2);
    p_params_->accelerometerCovariance = I3 * pow(settings_.accel_noise, 2);
    p_params_->integrationCovariance = I3 * 1e-8;
    p_params_->biasAccCovariance = I3 * pow(settings_.accel_bias_noise, 2);
    p_params_->biasOmegaCovariance = I3 * pow(settings_.gyro_bias_noise, 2);

    ISAM2Params isam_params;
    isam_params.relinearizeThreshold = 0.1;
    isam_params.relinearizeSkip = 50;  
    isam_params.factorization = ISAM2Params::CHOLESKY;
    isam_ = std::make_unique<ISAM2>(isam_params);

    prev_bias_ = imuBias::ConstantBias(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
    pim_ = std::make_unique<PreintegratedCombinedMeasurements>(p_params_, prev_bias_);
}

void FGOEngine::initialize(double t0, const Eigen::Vector3d& att, const Eigen::Vector3d& vn, const Eigen::Vector3d& pos,
                           const Eigen::Vector3d& bg_init, const Eigen::Vector3d& ba_init) {
    Rot3 R = Rot3::Ypr(att(2), att(0), att(1)); 
    Point3 P = pos; 
    prev_pose_ = Pose3(R, P);
    prev_vel_  = vn;
    prev_bias_ = imuBias::ConstantBias(ba_init, bg_init);
    
    pim_->resetIntegrationAndSetBias(prev_bias_);
    current_node_time_ = t0;

    // 初始先验
    auto noise_pose = noiseModel::Diagonal::Sigmas((gtsam::Vector(6) << 0.002, 0.002, 0.01, 0.1, 0.1, 0.1).finished()); 
    auto noise_vel  = noiseModel::Diagonal::Sigmas(Eigen::Vector3d(0.01, 0.01, 0.01));
    auto noise_bias = noiseModel::Diagonal::Sigmas((gtsam::Vector(6) << 
        1e-4, 1e-4, 1e-4,   // Acc Bias
        1e-5, 1e-5, 1e-5    // Gyro Bias
    ).finished());

    new_factors_.add(PriorFactor<Pose3>(X(0), prev_pose_, noise_pose));
    new_factors_.add(PriorFactor<Eigen::Vector3d>(V(0), prev_vel_, noise_vel));
    new_factors_.add(PriorFactor<imuBias::ConstantBias>(B(0), prev_bias_, noise_bias));

    new_values_.insert(X(0), prev_pose_);
    new_values_.insert(V(0), prev_vel_);
    new_values_.insert(B(0), prev_bias_);
    
    try {
        isam_->update(new_factors_, new_values_);
    } catch(const std::exception& e) {
        cerr << "❌ [Init Error] " << e.what() << endl;
    }
    new_factors_.resize(0); new_values_.clear();
    initialized_ = true;
    history_nodes_[t0] = 0;
    
    cout << "✅ [FGO] Initialized at t=" << t0 << endl;
    cout << "   [System] Vertical Damping ENABLED (Matching SINS behavior)" << endl;
}

bool FGOEngine::process_imu(const Eigen::Vector3d& gyro, const Eigen::Vector3d& acc, double dt, bool force_align, bool use_cai) {
    if (!initialized_) return false;
    if (gyro.hasNaN() || acc.hasNaN()) return false;
    
    pim_->integrateMeasurement(acc, gyro, dt);
    current_node_time_ += dt;

    if (pim_->deltaTij() >= settings_.node_interval - 1e-5) {
        key_idx_++;
        
        NavState prop = pim_->predict(NavState(prev_pose_, prev_vel_), prev_bias_);
        
        new_factors_.add(CombinedImuFactor(X(key_idx_-1), V(key_idx_-1), X(key_idx_), V(key_idx_), 
                                           B(key_idx_-1), B(key_idx_), *pim_));
        new_values_.insert(X(key_idx_), prop.pose());
        new_values_.insert(V(key_idx_), prop.velocity());
        new_values_.insert(B(key_idx_), prev_bias_);

        // -----------------------------------------------------------------------
        // [新增] 高度阻尼 (Vertical Damping) 逻辑
        // -----------------------------------------------------------------------
        // 即使 force_align=false (不锁水平速度)，我们也必须锁垂直速度，
        // 以对齐 SINS 的 vertical_damping_mode = 1。
        
        // 1. 垂直速度约束 (Vz -> 0)
        // Sigma配置: 水平=1000(完全不约束), 垂直=0.01(强约束)
        auto noise_vel_damping = noiseModel::Diagonal::Sigmas((gtsam::Vector(3) << 1000.0, 1000.0, 0.01).finished());
        new_factors_.add(PriorFactor<Eigen::Vector3d>(V(key_idx_), Eigen::Vector3d::Zero(), noise_vel_damping));

        // 2. 高度约束 (Pz -> 0, assuming local frame)
        // Sigma配置: 旋转=1000, 水平位置=1000, 垂直位置=0.01
        // 注意: GTSAM Pose3 sigmas 是 [Rot(3), Pos(3)] 即 [rx,ry,rz, x,y,z]
        auto noise_height_damping = noiseModel::Diagonal::Sigmas((gtsam::Vector(6) << 
            1000.0, 1000.0, 1000.0,  // Rotation (Unconstrained)
            1000.0, 1000.0, 0.01     // Position (Z locked)
        ).finished());
        
        // 这里的 Pose3() 仅用于提供 Z=0 的基准，X/Y/Rot 不会被约束因为 Sigma 很大
        new_factors_.add(PriorFactor<Pose3>(X(key_idx_), Pose3(), noise_height_damping));

        // -----------------------------------------------------------------------
        // [原有] 全向零速修正 (ZUPT)
        // -----------------------------------------------------------------------
        if (force_align) {
            auto noise_vel = noiseModel::Isotropic::Sigma(3, 0.001); 
            auto noise_pos = noiseModel::Diagonal::Sigmas((gtsam::Vector(6) << 100,100,100, 0.05,0.05,0.05).finished());
            
            new_factors_.add(PriorFactor<Eigen::Vector3d>(V(key_idx_), Eigen::Vector3d::Zero(), noise_vel));
            new_factors_.add(PriorFactor<Pose3>(X(key_idx_), Pose3(prev_pose_.rotation(), Point3(0,0,0)), noise_pos));
        } 
        
        try {
            isam_->update(new_factors_, new_values_);
            Values result = isam_->calculateEstimate();
            prev_pose_ = result.at<Pose3>(X(key_idx_));
            prev_vel_  = result.at<Eigen::Vector3d>(V(key_idx_));
            prev_bias_ = result.at<imuBias::ConstantBias>(B(key_idx_));
            
            history_nodes_[current_node_time_] = key_idx_;
            
        } catch(const gtsam::IndeterminantLinearSystemException& e) {
            cerr << "\n🔥 [CRASH] " << e.what() << endl;
            return false;
        } catch(const std::exception& e) {
            cerr << "\n🔥 [CRASH] " << e.what() << endl;
            return false;
        }

        new_factors_.resize(0);
        new_values_.clear();
        pim_->resetIntegration();
        return true;
    }
    return false;
}

void FGOEngine::add_cai_measurement(double t_end, const Eigen::Vector3d& w_avg) {
    auto it_curr = history_nodes_.lower_bound(t_end - 0.01);
    if (it_curr == history_nodes_.end()) return;
    
    double t_start = t_end - 2.0; 
    auto it_prev = history_nodes_.lower_bound(t_start - 0.01);
    if (it_prev == history_nodes_.end()) return;

    size_t key_curr = it_curr->second;
    size_t key_prev = it_prev->second;
    
    if (key_curr <= key_prev) return;

    double dt = t_end - t_start; 
    Eigen::Vector3d dtheta = w_avg * dt; 
    Rot3 dR = Rot3::Expmap(dtheta);

    Pose3 between_pose(dR, Point3(0,0,0));

    gtsam::Vector6 sigmas;
    sigmas << settings_.cai_rot_noise, settings_.cai_rot_noise, settings_.cai_rot_noise, 
              1000.0, 1000.0, 1000.0;
    
    auto noise = noiseModel::Diagonal::Sigmas(sigmas);
    
    new_factors_.add(BetweenFactor<Pose3>(X(key_prev), X(key_curr), between_pose, noise));
}

NavStateResult FGOEngine::get_result() {
    return {current_node_time_, prev_pose_, prev_vel_, prev_bias_};
}