#include "cai_sim.h"
#include <iostream>
#include <cmath>

namespace cai {

using namespace Eigen;

AtomicGyroSimulator::AtomicGyroSimulator(const Vector3d& pos_rad, const GLV& glv,
                                         const CAIParams& params)
    : params_(params), glv_(glv), rng_(params.seed), noise_dist_(0.0, 1.0)
{
    // 计算地球自转在 n 系的投影
    double lat = pos_rad(0);
    wie_n_ = Vector3d(0, glv_.wie * cos(lat), glv_.wie * sin(lat));
    wie_b_.setZero();
    
    // 计算噪声参数
    // ARW (deg/√h) -> 角度噪声 (rad/cycle)
    double t_cycle_hours = params_.T_cycle / 3600.0;
    double angle_noise_deg = params_.arw_dpsh * sqrt(t_cycle_hours);
    angle_noise_rad_ = angle_noise_deg * glv_.deg;
    
    // 计算零偏 (deg/h -> rad/s)
    double bias_rad_s = (params_.bias_dph * glv_.deg) / 3600.0;
    bias_vec_ = Vector3d(bias_rad_s, bias_rad_s, bias_rad_s);
}

void AtomicGyroSimulator::Init(const Vector3d& att_rad) {
    Matrix3d Cnb = INSMath::a2mat(att_rad);
    Init(Cnb);
}

void AtomicGyroSimulator::Init(const Matrix3d& Cnb) {
    wie_b_ = Cnb.transpose() * wie_n_;
    initialized_ = true;
}

void AtomicGyroSimulator::Init(const Quaterniond& qnb) {
    Matrix3d Cnb = INSMath::q2mat(qnb);
    Init(Cnb);
}

void AtomicGyroSimulator::SetWieB(const Vector3d& wie_b) {
    wie_b_ = wie_b;
    initialized_ = true;
}

Vector3d AtomicGyroSimulator::Measure() {
    if (!initialized_) {
        std::cerr << "[CAI] Warning: Not initialized, returning zero" << std::endl;
        return Vector3d::Zero();
    }
    return wie_b_ + bias_vec_;
}

Vector3d AtomicGyroSimulator::Measure(const Matrix3d& Cnb) {
    Vector3d wie_b = Cnb.transpose() * wie_n_;
    return wie_b + bias_vec_;
}

Vector3d AtomicGyroSimulator::Measure(const Quaterniond& qnb) {
    Matrix3d Cnb = INSMath::q2mat(qnb);
    return Measure(Cnb);
}

Vector3d AtomicGyroSimulator::GetAngleNoise() {
    return Vector3d(noise_dist_(rng_), noise_dist_(rng_), noise_dist_(rng_)) * angle_noise_rad_;
}

Vector3d AtomicGyroSimulator::GetTrueWieB(const Matrix3d& Cnb) const {
    return Cnb.transpose() * wie_n_;
}

Vector3d AtomicGyroSimulator::GetTrueWieB(const Vector3d& att_rad) const {
    Matrix3d Cnb = INSMath::a2mat(att_rad);
    return Cnb.transpose() * wie_n_;
}

void AtomicGyroSimulator::PrintConfig() const {
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "[Atomic Gyro Config]" << std::endl;
    std::cout << "  T_cycle:    " << params_.T_cycle << " s" << std::endl;
    std::cout << "  ARW:        " << params_.arw_dpsh << " deg/√h" << std::endl;
    std::cout << "  Bias:       " << params_.bias_dph << " deg/h" << std::endl;
    std::cout << "  Angle noise/cycle: " << angle_noise_rad_ / glv_.deg * 3600 << " arcsec" << std::endl;
    std::cout << "  wie_n (deg/h): " << (wie_n_ / glv_.deg * 3600).transpose() << std::endl;
    if (initialized_) {
        std::cout << "  wie_b (deg/h): " << (wie_b_ / glv_.deg * 3600).transpose() << std::endl;
    }
    std::cout << "---------------------------------------" << std::endl;
}

} // namespace cai
