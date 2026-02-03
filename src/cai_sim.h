#ifndef CAI_SIM_H
#define CAI_SIM_H

/**
 * @file cai_sim.h
 * @brief 原子陀螺(CAI/Atomic Gyro)仿真器 - 统一接口
 * 
 * 物理模型：
 *   - 测量地球自转在体坐标系的投影 (wie_b)
 *   - 包含散粒噪声 (ARW) 和零偏不稳定性
 *   - 干涉周期 T_cycle = 2.0s
 * 
 * 参数说明：
 *   - ARW = 2e-4 deg/√h (角度随机游走)
 *   - Bias Instability = 1e-5 deg/h (零偏不稳定性)
 */

#include <Eigen/Dense>
#include <random>
#include "glv.h"
#include "ins_math.h"

namespace cai {

// 原子陀螺物理参数
struct CAIParams {
    double T_cycle = 2.0;           // 干涉周期 (s)
    double arw_dpsh = 2.0e-5;       // ARW (deg/√h)
    double bias_dph = 1.0e-5;       // 零偏不稳定性 (deg/h)
    double bias_ug = 0.05;   // 零偏稳定性 (50ng = 0.05ug)
    double vrw_ug = 10.0;    // 速度随机游走/白噪声 (ug)
    int seed = 123;          // 随机种子
};

/**
 * @brief 原子陀螺仿真器
 * 
 * 使用方式：
 *   1. 构造时传入位置和参数
 *   2. 调用 Init() 设置初始姿态（用于计算 wie_b）
 *   3. 调用 Measure() 获取测量值
 */
class AtomicGyroSimulator {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    // 构造函数
    AtomicGyroSimulator(const Eigen::Vector3d& pos_rad, const GLV& glv, 
                        const CAIParams& params = CAIParams());
    
    // 用姿态初始化（计算固定的 wie_b）
    void Init(const Eigen::Vector3d& att_rad);
    void Init(const Eigen::Matrix3d& Cnb);
    void Init(const Eigen::Quaterniond& qnb);
    
    // 直接设置 wie_b（用于外部计算好的情况）
    void SetWieB(const Eigen::Vector3d& wie_b);
    
    // 获取测量值（使用初始化时的固定 wie_b）
    Eigen::Vector3d Measure();
    
    // 获取测量值（使用实时姿态计算 wie_b）
    Eigen::Vector3d Measure(const Eigen::Matrix3d& Cnb);
    Eigen::Vector3d Measure(const Eigen::Quaterniond& qnb);
    
    // 获取角度噪声（每周期）
    Eigen::Vector3d GetAngleNoise();
    
    // 获取 wie_b 真值
    Eigen::Vector3d GetTrueWieB() const { return wie_b_; }
    Eigen::Vector3d GetTrueWieB(const Eigen::Matrix3d& Cnb) const;
    Eigen::Vector3d GetTrueWieB(const Eigen::Vector3d& att_rad) const;
    
    // 获取 wie_n
    Eigen::Vector3d GetWieN() const { return wie_n_; }
    
    // 获取参数
    const CAIParams& GetParams() const { return params_; }
    double GetAngleNoiseStd() const { return angle_noise_rad_; }
    
    // 打印配置信息
    void PrintConfig() const;

private:
    CAIParams params_;
    GLV glv_;
    
    Eigen::Vector3d wie_n_;          // 地球自转在 n 系 (rad/s)
    Eigen::Vector3d wie_b_;          // 地球自转在 b 系 (rad/s) - 初始化时固定
    Eigen::Vector3d bias_vec_;       // 零偏向量 (rad/s)
    double angle_noise_rad_;         // 角度噪声标准差 (rad/cycle)
    
    std::mt19937 rng_;
    std::normal_distribution<double> noise_dist_;
    
    bool initialized_ = false;
};

/**
 * @brief 原子加速度计仿真器
 * * 物理模型：
 * - 测量比力 (Specific Force)
 * - 使用对准锁定的姿态计算理想参考比力 fb_ref = Cnb^T * (-g^n)
 * - 误差模型：Bias (ug) + White Noise (ug)
 */
class AtomicAccSimulator {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // 构造函数
    AtomicAccSimulator(const Eigen::Vector3d& pos_rad, const GLV& glv, 
                       const CAIParams& params = CAIParams());

    // 初始化：传入对准锁定的姿态 (Attitude Locked)
    // 这一步非常关键：它计算出所谓的“完美比力”，用于抵消 FOG 加计的零偏
    void Init(const Eigen::Vector3d& att_locked);

    // 获取测量值 (m/s^2)
    Eigen::Vector3d Measure();

private:
    CAIParams params_;
    GLV glv_;
    Eigen::Vector3d fn_ref_;     // n系参考比力 (指向天向)
    Eigen::Vector3d fb_ref_;     // b系参考比力 (初始化后固定)
    Eigen::Vector3d bias_vec_;   // 仿真零偏
    
    std::mt19937 rng_;
    std::normal_distribution<double> noise_dist_;
    
    bool initialized_ = false;
};


} // namespace cai

#endif // CAI_SIM_H
