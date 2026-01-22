#ifndef CAI_SIM_H
#define CAI_SIM_H

#include <Eigen/Dense>
#include <random>
#include <iostream>
#include "glv.h"

using namespace Eigen;
using namespace std;

class CAIGSimulator {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    CAIGSimulator();
    
    void Init(double start_time_sec);

    bool Update(double t_curr, const Eigen::Quaterniond& qnb_curr, Eigen::Vector3d& out_omega);
    
private:
    // --- 内部状态 ---
    double start_time_;       
    double cycle_ref_time_;   
    bool is_running_;         
    bool is_integrating_;     

    // --- 缓存数据 ---
    Eigen::Quaterniond qnb_start_latch_; 
    
    // --- 噪声与参数对象 ---
    std::default_random_engine generator_;
    std::normal_distribution<double> noise_dist_;
    Eigen::Vector3d bias_vec_;
    GLV glv_; // 确保 GLV 有默认构造函数初始化常数

    // --- 物理常量 ---
    const double T_CYCLE  = 2.0; 
    const double T_DEAD   = 0.4; 
    const double T_INTERF = 1.6; 
};

#endif // CAI_SIM_H