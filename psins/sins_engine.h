#ifndef SINS_ENGINE_H
#define SINS_ENGINE_H

#include <vector>
#include <iostream>
#include <fstream>
#include "glv.h"
#include "earth.h"
#include "ins_state.h"
#include "ins_wrapper.h" // 包含 RunAlignment 和 AlignConfig

// 引擎工作状态
enum class EngineState {
    IDLE,       // 未初始化
    ALIGNING,   // 对准中
    NAVIGATING, // 导航中
    FAULT       // 故障/炸机
};

class SinsEngine {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // --- 核心成员 ---
    GLV glv;
    Earth eth;
    INSState ins;       // 纯惯导状态
    AlignmentEngine aligner; // 对准引擎 (新增)
    
    EngineState state;
    double align_duration;   // 对准需持续的时间 (s)
    AlignConfig align_cfg;   // 对准参数配置

    // --- 构造与配置 ---
    SinsEngine(double ts = 1.0/400.0);

    /**
     * @brief 初始化引擎配置
     * @param cfg 对准参数（包含参考位置、参考姿态、噪声参数等）
     * @param align_time_sec 对准持续时间 (秒)
     */
    void Init(const AlignConfig& cfg, double align_time_sec);

    // --- 核心驱动 ---
    /**
     * @brief 智能单步推演：输入一帧数据，内部自动处理对准/导航切换
     * @param wm 陀螺角增量 (rad)
     * @param vm 加计速度增量 (m/s)
     * @param t  全局时间戳 (s)
     */
    void Step(const Eigen::Vector3d& wm, const Eigen::Vector3d& vm, double t);

    // --- 输出接口 ---
    Eigen::Vector3d GetAttDeg() const; // 欧拉角 (deg)
    Eigen::Vector3d GetVel() const;    // 速度 (m/s)
    Eigen::Vector3d GetPosDeg() const; // 位置 (deg, deg, m)
    Eigen::Vector3d GetBiasGyro() const; // 陀螺零偏 (rad/s)
    Eigen::Vector3d GetBiasAcc() const;  // 加计零偏 (m/s^2)
    Eigen::Quaterniond GetQnb() const;
    // --- 外部干预 ---
    void InjectBias(const Eigen::Vector3d& db_gyro, const Eigen::Vector3d& db_acc);
};

#endif // SINS_ENGINE_H