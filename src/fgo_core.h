#pragma once
#include <gtsam/navigation/CombinedImuFactor.h>
#include <gtsam/nonlinear/ISAM2.h> // ✅ 使用最稳的 ISAM2
#include <gtsam/nonlinear/Values.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/slam/PriorFactor.h>
#include <gtsam/slam/BetweenFactor.h>
#include <boost/shared_ptr.hpp> 
#include <Eigen/Dense> 

#include "../psins/glv.h"

using namespace gtsam;

struct FGOSettings {
    GLV glv_;
    double fog_dt = 1.0 / 400.0;
    double node_interval = 1.0; 
    double gravity = glv_.g;

    double gyro_noise = 1e-6;      
    double accel_noise = 1e-4;     
    double gyro_bias_noise = 1e-8; 
    double accel_bias_noise = 1e-6;
    
    double cai_rot_noise = 1e-5; 
};

struct NavStateResult {
    double time;
    Pose3 pose;
    Eigen::Vector3d vel; 
    imuBias::ConstantBias bias;
};

class FGOEngine {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW 

    FGOEngine(const FGOSettings& settings);
    
    void initialize(double t0, const Eigen::Vector3d& att, const Eigen::Vector3d& vn, const Eigen::Vector3d& pos,
                    const Eigen::Vector3d& bg_init, const Eigen::Vector3d& ba_init);

    bool process_psins_delta(const Pose3& dPose_psins, const Eigen::Vector3d& vel_curr_psins, double t_curr);
    void add_cai_measurement(double t_end, const Eigen::Vector3d& w_avg);
    NavStateResult get_result();

private:
    FGOSettings settings_;
    std::unique_ptr<ISAM2> isam_; // ✅ 稳如老狗的引擎
    
    Pose3 prev_pose_;
    Eigen::Vector3d prev_vel_; 
    imuBias::ConstantBias prev_bias_;

    NonlinearFactorGraph new_factors_;
    Values new_values_;
    
    uint64_t key_idx_ = 0;
    double current_node_time_ = 0.0;
    bool initialized_ = false;

    std::map<double, uint64_t> history_nodes_; 
};