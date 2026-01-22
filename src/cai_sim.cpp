#include "cai_sim.h"
#include <iostream> 

using namespace Eigen;

CAIGSimulator::CAIGSimulator() 
    : start_time_(-1.0), cycle_ref_time_(-1.0), 
      is_running_(false), is_integrating_(false) {
    // (1) 角度随机游走 (ARW) / 散粒噪声
    // 设定目标：2e-4 deg/sqrt(h) -> 这是一个高性能战术级/导航级的水平
    double arw_dpsh = 2.0e-4; 

    // (2) 零偏不稳定性 (Bias Instability)
    // 原子陀螺的强项：长期稳定性极好，通常比 FOG 好 1-2 个数量级
    // 设定目标：5e-5 deg/h
    double bias_instability_dph = 5.0e-5;       

    // ==========================================
    // 2. 数值转换 (Convert to SI units for Simulation)
    // ==========================================

    // --- 计算单次采样噪声标准差 (Sigma) ---
    // 公式: sigma_rate = ARW / sqrt(T_cycle)
    // 单位换算:
    //   T_cycle = 2.0 秒 = 2.0 / 3600 小时
    //   sigma (deg/h) = 2e-4 / sqrt(2/3600) 
    //                 = 2e-4 / 0.02357 ≈ 0.00848 deg/h
    double t_cycle_hours = T_CYCLE / 3600.0;
    double sigma_deg_h = arw_dpsh / sqrt(t_cycle_hours);
    
    // 最终转为 rad/s 供 update 函数使用
    double sigma_rad_s = (sigma_deg_h * glv_.deg) / 3600.0;

    // --- 计算零偏 ---
    double bias_rad_s = (bias_instability_dph * glv_.deg) / 3600.0;
    
    // ==========================================
    // 3. 初始化生成器
    // ==========================================
    noise_dist_ = std::normal_distribution<double>(0.0, sigma_rad_s);
    bias_vec_ << bias_rad_s, bias_rad_s, bias_rad_s;
    
    // 打印参数确认 (调试用)
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "[CAIG Physics Config]" << std::endl;
    std::cout << "  Target ARW:   " << arw_dpsh << " deg/sqrt(h)" << std::endl;
    std::cout << "  Cycle Time:   " << T_CYCLE << " s" << std::endl;
    std::cout << "  Calc Sigma:   " << sigma_deg_h << " deg/h (Per Sample)" << std::endl;
    std::cout << "  Calc Bias:    " << bias_instability_dph << " deg/h" << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    qnb_start_latch_.setIdentity();
}

void CAIGSimulator::Init(double start_time_sec) {
    start_time_ = start_time_sec;
    cycle_ref_time_ = start_time_sec; 
    is_running_ = false;
    is_integrating_ = false;
    qnb_start_latch_.setIdentity();
    
    std::cout << "[CAIG] Initialized. Start t = " << start_time_ << " s." << std::endl;
}

bool CAIGSimulator::Update(double t_curr, const Quaterniond& qnb_curr, Vector3d& out_omega) {
    if (t_curr < start_time_) return false;

    if (!is_running_) {
        is_running_ = true;
        cycle_ref_time_ = start_time_; 
    }

    double t_local = t_curr - cycle_ref_time_;

    // [阶段 A]: 死区 (0 ~ 0.4s)
    if (t_local < T_DEAD) {
        is_integrating_ = false; 
        return false;
    }

    // [阶段 B]: 干涉积分区 (0.4s ~ 2.0s)
    else if (t_local < T_CYCLE) {
        if (!is_integrating_) {
            qnb_start_latch_ = qnb_curr;
            is_integrating_ = true;
        }
        return false;
    }

    // [阶段 C]: 周期结束
    else {
        // q_diff = q_start^{-1} * q_end (Body frame increment)
        Quaterniond q_diff = qnb_start_latch_.conjugate() * qnb_curr;
        q_diff.normalize(); // [Fix] 增加归一化，防止数值误差积累

        Eigen::AngleAxisd rotation_vector(q_diff);
        double angle = rotation_vector.angle();
        
        // 规范化角度到 -pi ~ pi
        if (angle > M_PI) angle -= 2.0 * M_PI;
        else if (angle < -M_PI) angle += 2.0 * M_PI;

        Vector3d w_avg;
        // [Fix] 增加小角度保护
        if (std::abs(angle) < 1e-9) {
            w_avg.setZero();
        } else {
            w_avg = (rotation_vector.axis() * angle) / T_INTERF; 
        }

        double n1 = noise_dist_(generator_);
        double n2 = noise_dist_(generator_);
        double n3 = noise_dist_(generator_);
        Vector3d noise(n1, n2, n3);

        out_omega = w_avg + bias_vec_ + noise;

        cycle_ref_time_ += T_CYCLE; 
        is_integrating_ = false; 

        return true; 
    }
}