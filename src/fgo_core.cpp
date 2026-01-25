#include "fgo_core.h"
#include <iostream>
#include <boost/make_shared.hpp>
#include <gtsam/linear/linearExceptions.h> 

using namespace std;
using namespace gtsam;
using namespace gtsam::symbol_shorthand;

FGOEngine::FGOEngine(const FGOSettings& settings) : settings_(settings) {
    ISAM2Params isam_params;
    isam_params.relinearizeThreshold = 0.1;
    isam_params.relinearizeSkip = 10; 
    isam_params.factorization = ISAM2Params::CHOLESKY;
    isam_ = std::make_unique<ISAM2>(isam_params);

    prev_bias_ = imuBias::ConstantBias(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
}

void FGOEngine::initialize(double t0, const Eigen::Vector3d& att, const Eigen::Vector3d& vn, const Eigen::Vector3d& pos,
                           const Eigen::Vector3d& bg_init, const Eigen::Vector3d& ba_init) {
    Rot3 R = Rot3::Ypr(att(2), att(0), att(1)); 
    Point3 P = pos; 
    prev_pose_ = Pose3(R, P);
    prev_vel_  = vn;
    prev_bias_ = imuBias::ConstantBias(ba_init, bg_init);
    
    current_node_time_ = t0;

    auto noise_pose = noiseModel::Diagonal::Sigmas((gtsam::Vector(6) << 0.002, 0.002, 0.01, 0.1, 0.1, 0.1).finished()); 
    auto noise_vel  = noiseModel::Diagonal::Sigmas(Eigen::Vector3d(0.01, 0.01, 0.01));
    auto noise_bias = noiseModel::Diagonal::Sigmas((gtsam::Vector(6) << 1e-4, 1e-4, 1e-4, 1e-5, 1e-5, 1e-5).finished());

    new_factors_.add(PriorFactor<Pose3>(X(0), prev_pose_, noise_pose));
    new_factors_.add(PriorFactor<Eigen::Vector3d>(V(0), prev_vel_, noise_vel));
    new_factors_.add(PriorFactor<imuBias::ConstantBias>(B(0), prev_bias_, noise_bias));

    new_values_.insert(X(0), prev_pose_);
    new_values_.insert(V(0), prev_vel_);
    new_values_.insert(B(0), prev_bias_);
    
    try {
        isam_->update(new_factors_, new_values_);
    } catch(...) {}
    new_factors_.resize(0); new_values_.clear();
    initialized_ = true;
    history_nodes_[t0] = 0;
    
    cout << "✅ [FGO] Initialized (Engine: Stable ISAM2)" << endl;
}

bool FGOEngine::process_psins_delta(const Pose3& dPose_psins, const Eigen::Vector3d& vel_curr_psins, double t_curr) {
    if (!initialized_) return false;
    
    current_node_time_ = t_curr;
    key_idx_++;
    
    auto noise_pose = noiseModel::Diagonal::Sigmas((gtsam::Vector(6) << 1e-5, 1e-5, 1e-5, 1e-2, 1e-2, 1e-2).finished());
    new_factors_.add(BetweenFactor<Pose3>(X(key_idx_-1), X(key_idx_), dPose_psins, noise_pose));
    
    auto noise_vel = noiseModel::Diagonal::Sigmas(Eigen::Vector3d(1e-3, 1e-3, 1e-3));
    new_factors_.add(PriorFactor<Eigen::Vector3d>(V(key_idx_), vel_curr_psins, noise_vel));
    
    auto noise_bias = noiseModel::Diagonal::Sigmas((gtsam::Vector(6) << 1e-4,1e-4,1e-4, 1e-5,1e-5,1e-5).finished());
    new_factors_.add(BetweenFactor<imuBias::ConstantBias>(B(key_idx_-1), B(key_idx_), imuBias::ConstantBias(), noise_bias));
    
    Pose3 pred_pose = prev_pose_.compose(dPose_psins);
    new_values_.insert(X(key_idx_), pred_pose);
    new_values_.insert(V(key_idx_), vel_curr_psins);
    new_values_.insert(B(key_idx_), prev_bias_);
    
    try {
        isam_->update(new_factors_, new_values_);
        
        Values result = isam_->calculateEstimate();
        prev_pose_ = result.at<Pose3>(X(key_idx_));
        prev_vel_  = result.at<Eigen::Vector3d>(V(key_idx_));
        prev_bias_ = result.at<imuBias::ConstantBias>(B(key_idx_));
        
        history_nodes_[current_node_time_] = key_idx_;
        
        // 清理旧索引
        auto it = history_nodes_.begin();
        while(it != history_nodes_.end() && it->first < (t_curr - 10.0)) {
            it = history_nodes_.erase(it);
        }

    } catch(...) {
        cerr << "❌ [Update Error]" << endl;
        return false;
    }

    new_factors_.resize(0);
    new_values_.clear();
    return true;
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
              1.0e10, 1.0e10, 1.0e10; 
    auto noise = noiseModel::Diagonal::Sigmas(sigmas);
    
    new_factors_.add(BetweenFactor<Pose3>(X(key_prev), X(key_curr), between_pose, noise));
}

NavStateResult FGOEngine::get_result() {
    return {current_node_time_, prev_pose_, prev_vel_, prev_bias_};
}


// #include "fgo_core.h"
// #include <iostream>
// #include <boost/make_shared.hpp>
// #include <gtsam/linear/linearExceptions.h> 

// using namespace std;
// using namespace gtsam;
// using namespace gtsam::symbol_shorthand;

// FGOEngine::FGOEngine(const FGOSettings& settings) : settings_(settings) {
//     ISAM2Params isam_params;
//     isam_params.relinearizeThreshold = 0.1;
//     // [关键优化] 每10次更新才重组一次矩阵，极大提速！
//     // 这让我们能在 100Hz 下跑 ISAM2 而不卡顿
//     isam_params.relinearizeSkip = 10; 
//     isam_params.factorization = ISAM2Params::CHOLESKY;
//     isam_ = std::make_unique<ISAM2>(isam_params);

//     prev_bias_ = imuBias::ConstantBias(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
// }

// void FGOEngine::initialize(double t0, const Eigen::Vector3d& att, const Eigen::Vector3d& vn, const Eigen::Vector3d& pos,
//                            const Eigen::Vector3d& bg_init, const Eigen::Vector3d& ba_init) {
//     Rot3 R = Rot3::Ypr(att(2), att(0), att(1)); 
//     Point3 P = pos; 
//     prev_pose_ = Pose3(R, P);
//     prev_vel_  = vn;
//     prev_bias_ = imuBias::ConstantBias(ba_init, bg_init);
    
//     current_node_time_ = t0;

//     auto noise_pose = noiseModel::Diagonal::Sigmas((gtsam::Vector(6) << 0.002, 0.002, 0.01, 0.1, 0.1, 0.1).finished()); 
//     auto noise_vel  = noiseModel::Diagonal::Sigmas(Eigen::Vector3d(0.01, 0.01, 0.01));
//     auto noise_bias = noiseModel::Diagonal::Sigmas((gtsam::Vector(6) << 1e-4, 1e-4, 1e-4, 1e-5, 1e-5, 1e-5).finished());

//     new_factors_.add(PriorFactor<Pose3>(X(0), prev_pose_, noise_pose));
//     new_factors_.add(PriorFactor<Eigen::Vector3d>(V(0), prev_vel_, noise_vel));
//     new_factors_.add(PriorFactor<imuBias::ConstantBias>(B(0), prev_bias_, noise_bias));

//     new_values_.insert(X(0), prev_pose_);
//     new_values_.insert(V(0), prev_vel_);
//     new_values_.insert(B(0), prev_bias_);
    
//     try {
//         isam_->update(new_factors_, new_values_);
//     } catch(...) {}
//     new_factors_.resize(0); new_values_.clear();
//     initialized_ = true;
//     history_nodes_[t0] = 0;
    
//     cout << "✅ [FGO] Initialized (Engine: Stable ISAM2)" << endl;
// }

// bool FGOEngine::process_psins_delta(const Pose3& dPose_psins, const Eigen::Vector3d& vel_curr_psins, double t_curr) {
//     if (!initialized_) return false;
    
//     current_node_time_ = t_curr;
//     key_idx_++;
    
//     auto noise_pose = noiseModel::Diagonal::Sigmas((gtsam::Vector(6) << 1e-5, 1e-5, 1e-5, 1e-2, 1e-2, 1e-2).finished());
//     new_factors_.add(BetweenFactor<Pose3>(X(key_idx_-1), X(key_idx_), dPose_psins, noise_pose));
    
//     auto noise_vel = noiseModel::Diagonal::Sigmas(Eigen::Vector3d(1e-3, 1e-3, 1e-3));
//     new_factors_.add(PriorFactor<Eigen::Vector3d>(V(key_idx_), vel_curr_psins, noise_vel));
    
//     auto noise_bias = noiseModel::Diagonal::Sigmas((gtsam::Vector(6) << 1e-4,1e-4,1e-4, 1e-5,1e-5,1e-5).finished());
//     new_factors_.add(BetweenFactor<imuBias::ConstantBias>(B(key_idx_-1), B(key_idx_), imuBias::ConstantBias(), noise_bias));
    
//     Pose3 pred_pose = prev_pose_.compose(dPose_psins);
//     new_values_.insert(X(key_idx_), pred_pose);
//     new_values_.insert(V(key_idx_), vel_curr_psins);
//     new_values_.insert(B(key_idx_), prev_bias_);
    
//     try {
//         // [核心] ISAM2 更新
//         isam_->update(new_factors_, new_values_);
        
//         Values result = isam_->calculateEstimate();
//         prev_pose_ = result.at<Pose3>(X(key_idx_));
//         prev_vel_  = result.at<Eigen::Vector3d>(V(key_idx_));
//         prev_bias_ = result.at<imuBias::ConstantBias>(B(key_idx_));
        
//         history_nodes_[current_node_time_] = key_idx_;
        
//         // [内存管理] 仅为了查找 CAI 关联，我们可以安全清理 map
//         // 这不是边缘化，只是清理我们的索引表，完全安全！
//         auto it = history_nodes_.begin();
//         while(it != history_nodes_.end() && it->first < (t_curr - 10.0)) {
//             it = history_nodes_.erase(it);
//         }

//     } catch(...) {
//         cerr << "❌ [Update Error]" << endl;
//         return false;
//     }

//     new_factors_.resize(0);
//     new_values_.clear();
//     return true;
// }

// void FGOEngine::add_cai_measurement(double t_end, const Eigen::Vector3d& w_avg) {
//     auto it_curr = history_nodes_.lower_bound(t_end - 0.01);
//     if (it_curr == history_nodes_.end()) return;
    
//     double t_start = t_end - 2.0; 
//     auto it_prev = history_nodes_.lower_bound(t_start - 0.01);
//     if (it_prev == history_nodes_.end()) return;

//     size_t key_curr = it_curr->second;
//     size_t key_prev = it_prev->second;
//     if (key_curr <= key_prev) return;

//     double dt = t_end - t_start; 
//     Eigen::Vector3d dtheta = w_avg * dt; 
//     Rot3 dR = Rot3::Expmap(dtheta);
//     Pose3 between_pose(dR, Point3(0,0,0));

//     gtsam::Vector6 sigmas;
//     sigmas << settings_.cai_rot_noise, settings_.cai_rot_noise, settings_.cai_rot_noise, 
//               1.0e10, 1.0e10, 1.0e10; 
//     auto noise = noiseModel::Diagonal::Sigmas(sigmas);
    
//     new_factors_.add(BetweenFactor<Pose3>(X(key_prev), X(key_curr), between_pose, noise));
// }

// NavStateResult FGOEngine::get_result() {
//     return {current_node_time_, prev_pose_, prev_vel_, prev_bias_};
// }