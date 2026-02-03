// #ifndef FGO_MODULE_H
// #define FGO_MODULE_H

// #include <Eigen/Dense>
// #include <vector>
// #include <iostream>
// #include "glv.h"
// #include "ins_math.h"

// using namespace Eigen;
// using namespace std;

// // ==========================================
// // [FGO] 因子图优化与预积分模块
// // 特性：实时预积分、雅可比维护、滑动窗口优化
// // ==========================================

// class FGO_Optimizer {
// public:
//     EIGEN_MAKE_ALIGNED_OPERATOR_NEW

//     struct FGOState {
//         Vector3d bg; // 陀螺零偏
//         Vector3d ba; // 加计零偏
//     };

//     // 预积分数据结构 (Pre-integrated Measurements)
//     // 保存了窗口内的相对运动量和对零偏的敏感度(Jacobian)
//     struct PreIntegration {
//         double sum_dt = 0;
//         Vector3d delta_theta = Vector3d::Zero(); // 角度增量
//         Vector3d delta_vel   = Vector3d::Zero(); // 速度增量
//         Vector3d delta_pos   = Vector3d::Zero(); // 位置增量
        
//         // 雅可比矩阵: J_ba_p (位置对加计零偏的偏导数) 等
//         Matrix3d J_ba_v = Matrix3d::Zero(); 
//         Matrix3d J_ba_p = Matrix3d::Zero();
//         Matrix3d J_bg_theta = Matrix3d::Zero();
        
//         // 简单的协方差累积 (用于加权)
//         double cov_val = 0; 

//         void Reset() {
//             sum_dt = 0;
//             delta_theta.setZero(); delta_vel.setZero(); delta_pos.setZero();
//             J_ba_v.setZero(); J_ba_p.setZero(); J_bg_theta.setZero();
//             cov_val = 0;
//         }
//     };

//     FGO_Optimizer(double ts) : ts_(ts) {
//         pim_.Reset();
//     }

//     void Reset(const Vector3d& init_bg, const Vector3d& init_ba) {
//         curr_bg_ = init_bg;
//         curr_ba_ = init_ba;
//         pim_.Reset();
//     }

//     // ---------------------------------------------------------
//     // [实时部分] 400Hz 调用：预积分更新
//     // 这一步完全非阻塞，只做简单的累加和雅可比传播
//     // ---------------------------------------------------------
//     void Propagate(const Vector3d& wm, const Vector3d& vm, const Matrix3d& Cnb_curr) {
//         double dt = ts_;
        
//         // 1. 扣除当前零偏估计 (First-order approx)
//         Vector3d unbias_gyro = wm - curr_bg_ * dt; // 这里的wm是角度增量
//         Vector3d unbias_acc  = vm - curr_ba_ * dt; // 这里的vm是速度增量(比力积分)

//         // 2. 状态预积分 (简化版: 忽略地球自转 Coriolis，只看本体运动)
//         // delta_p += delta_v * dt + 0.5 * R * acc * dt^2
//         pim_.delta_pos += pim_.delta_vel * dt + 0.5 * Cnb_curr * unbias_acc * dt; 
//         pim_.delta_vel += Cnb_curr * unbias_acc;
//         pim_.delta_theta += unbias_gyro; // 简化为直接累加，实际应为四元数乘法

//         // 3. 雅可比传播 (核心：记录状态对零偏的敏感度)
//         // d(vel) / d(ba) = - R * dt
//         // d(pos) / d(ba) = - 0.5 * R * dt^2
//         pim_.J_ba_v -= Cnb_curr * dt;
//         pim_.J_ba_p -= 0.5 * Cnb_curr * dt * dt;
//         pim_.J_bg_theta -= Matrix3d::Identity() * dt;

//         pim_.sum_dt += dt;
//         pim_.cov_val += 1.0; // 简化协方差传播
//     }

//     // ---------------------------------------------------------
//     // [优化部分] 0.5Hz 调用：原子数据触发因子图求解
//     // ---------------------------------------------------------
//     void Solve(const Vector3d& atom_gyro_meas, const Vector3d& atom_acc_meas, double gain) {
//         // 构建因子图残差模型
//         // 观测方程: Meas = True_Motion + Bias + Noise
//         // 预积分方程: PIM = True_Motion
//         // 这里的逻辑稍微做个变换：我们直接比较 (IMU均值 - 原子测量) 来估算 Bias 误差
        
//         // FOG 平均观测值 (基于预积分结果)
//         Vector3d fog_gyro_mean = pim_.delta_theta / pim_.sum_dt;
//         Vector3d fog_acc_mean  = pim_.delta_vel / pim_.sum_dt; // 注意这里是用速度增量除以时间得到平均加速度

//         // 残差 r = z - Hx
//         // 我们的目标是找到 delta_bias，使得观测残差最小
//         Vector3d r_gyro = fog_gyro_mean - atom_gyro_meas;
//         Vector3d r_acc  = fog_acc_mean  - atom_acc_meas;

//         // 一步高斯-牛顿更新 (One-step Gauss-Newton)
//         // 假设 H = I (单位阵), R 极小 (原子极准)
//         // 则 delta_x = K * r, 其中 K 由增益 gain 控制
        
//         Vector3d delta_bg = r_gyro * gain;
//         Vector3d delta_ba = r_acc  * gain;

//         // 更新状态
//         curr_bg_ += delta_bg; // 累加修正量
//         curr_ba_ += delta_ba;

//         // 存储本轮的修正量，用于修正历史轨迹
//         last_delta_bg_ = delta_bg;
//         last_delta_ba_ = delta_ba;
        
//         // 优化结束，重置预积分器
//         pim_.Reset();
//     }

//     // ---------------------------------------------------------
//     // [回溯修正] 线性化修正历史数据 (Linear Correction)
//     // 区别于"重算"：这里利用雅可比直接把 Bias 的变化量映射到 Pos/Vel 变化量
//     // Pos_new = Pos_old + J_ba * delta_ba
//     // ---------------------------------------------------------
//     void CorrectTrajectory(std::vector<IMUData>& buffer, 
//                            std::vector<Vector3d>& pos_hist, 
//                            std::vector<Vector3d>& vel_hist) {
        
//         // 简单的线性衰减修正：离当前时刻越近，受本次Bias突变影响越小（累积时间短）
//         // 但为了工程实现简单且效果显著，我们采用统一的梯度修正
//         // 在严谨的 FGO 中，这里会遍历 buffer，根据每个时刻的 Jacobian * delta_bias 进行修正
        
//         // 这里演示核心思想：不重新运行 Step_Nav，直接修数值
//         double dt = ts_;
//         Vector3d acc_correction = last_delta_ba_; // 加计零偏的变化量
        
//         // 对过去窗口内的每个点，应用 "由于零偏变化导致的积分误差"
//         // 误差传播公式: err_pos(t) = 0.5 * acc_err * t^2
//         for (size_t i = 0; i < pos_hist.size(); ++i) {
//             double t = (i + 1) * dt; // 积分时间
            
//             // 修正位置: P_new = P_old - (0.5 * delta_ba * t^2)
//             // 注意符号：如果零偏估计变大了，说明之前减多了，或者反之。
//             // 这里 Bias 是加在测量值上的， INS 计算是 (Meas - Bias)。
//             // 如果 Bias 增加 delta，则 (Meas - Bias) 减小 delta。
//             // 积分项减小 0.5 * delta * t^2。
            
//             Vector3d d_pos = -0.5 * acc_correction * t * t; 
//             Vector3d d_vel = -acc_correction * t;

//             // 甚至不需要旋转矩阵 Rnb，因为这本身就是在 n 系的等效误差
//             // (这是简化处理，严谨处理需要乘 Rnb，但短时间 Rnb 变化不大)
            
//             pos_hist[i] += d_pos; // 直接修改缓存中的历史数据
//             vel_hist[i] += d_vel;
//         }
//     }

//     Vector3d GetBg() const { return curr_bg_; }
//     Vector3d GetBa() const { return curr_ba_; }

// private:
//     double ts_;
//     PreIntegration pim_;
//     Vector3d curr_bg_ = Vector3d::Zero();
//     Vector3d curr_ba_ = Vector3d::Zero();
    
//     Vector3d last_delta_bg_ = Vector3d::Zero();
//     Vector3d last_delta_ba_ = Vector3d::Zero();
// };

// #endif // FGO_MODULE_H