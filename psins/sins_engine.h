#ifndef SINS_ENGINE_H
#define SINS_ENGINE_H

#include <vector>
#include <Eigen/Dense>
#include "glv.h"
#include "earth.h"
#include "ins_math.h"
#include "ins_state.h"
#include "kf_state.h"

struct IMUData {
    double t;
    Eigen::Vector3d wm; // 增量 rad
    Eigen::Vector3d vm; // 增量 m/s
};

struct AlignResult {
    bool valid = false;
    double align_time = 0.0;     // 对准持续时长
    Eigen::Vector3d att;         // 姿态 (rad)
    Eigen::Vector3d vel;         // 速度
    Eigen::Vector3d pos;         // 位置
    Eigen::Vector3d eb;          // 陀螺零偏
    Eigen::Vector3d db;          // 加计零偏
};

struct KFConfig {
    Eigen::Vector3d phi_init_err = {0.1*M_PI/180, 0.1*M_PI/180, 1.0*M_PI/180}; 
    Eigen::Vector3d wvn_err      = {0.001, 0.001, 0.001};         

    double eb_sigma = 0.2 * M_PI/180 / 3600.0;             
    double db_sigma = 100.0 * 1e-6 * 9.78;      // 100 ug
    // Qk: 角速度随机游走 (ARW, 对应Allan图 tau=1 处的值)
    double web_psd  = 0.02 * M_PI/180 / 60.0;              
    double wdb_psd  = 10.0 * 1e-6 * 9.78;       // 10 ug/sqrt(Hz)
};

// 混合对准配置
struct HybridAlignConfig {
    double t_coarse = 60.0;         // 粗对准时长 (s)
    double t_fine = 3600.0;         // 精对准时长 (s) 
    double eb_sigma_allan = 0.003;  // Allan 零偏稳定性 (deg/h)
    bool verbose = true;            // 是否输出过程信息
};

// 混合对准结果
struct HybridAlignResult {
    bool valid = false;
    Eigen::Vector3d att;     // 最终姿态 (rad)
    Eigen::Vector3d eb;      // 几何均值零偏 (rad/s)
    Eigen::Vector3d db;      // 加计零偏
    double align_time;       // 对准总时长
    
    // 中间结果 (调试用)
    Eigen::Vector3d att_coarse;  // 粗对准姿态
    Eigen::Vector3d att_kf;      // KF 精对准姿态
    Eigen::Vector3d eb_raw;      // 原始计算零偏
    Eigen::Vector3d eb_scale;    // 缩放因子
};

class SinsEngine {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // 基础模块
    GLV glv;
    Earth eth;
    INSState ins;       
    KFConfig kf_cfg;

    // 配置备份
    AlignResult res_init;   // 初始默认值 (Sins_Init 设定)
    AlignResult res_coarse; // 粗对准产出
    AlignResult res_fine;   // 精对准产出 

    Eigen::Vector3d coarse_sum_wm, coarse_sum_vm;
    double coarse_timer;       

    KFAlignVN kf;
    Eigen::Quaterniond kf_qnb;
    Eigen::Vector3d kf_vn, kf_pos;
    struct { Eigen::Vector3d wm, vm; } kf_last_imu;
    bool kf_first_step;    
    Eigen::Vector3d kf_base_eb, kf_base_db;
    double scale_ratio = 1.0;
public:   
    SinsEngine(double ts = 1.0/400.0);

    void Sins_Init(const Eigen::Vector3d& init_pos, // [Lat, Lon, H]
                   const Eigen::Vector3d& init_vel, // [Ve, Vn, Vu]
                   const Eigen::Vector3d& init_att, // [Pitch, Roll, Yaw]
                   const Eigen::Vector3d& init_eb,  // [Gx, Gy, Gz]
                   const Eigen::Vector3d& init_db); // [Ax, Ay, Az]
    void SetKFConfig(const KFConfig& cfg);

    void Run_Coarse_Phase(const std::vector<IMUData>& data_chunk);
    void Run_Fine_Phase(const std::vector<IMUData>& data_chunk);
    void Run_Nav_Phase(const std::vector<IMUData>& data_chunk);    

    Eigen::Vector3d GetAttDeg() const;    // 输出: [Pitch, Roll, Yaw] (deg)
    Eigen::Vector3d GetVel() const;       // 输出: [Ve, Vn, Vu] (m/s)
    Eigen::Vector3d GetPosDeg() const;    // 输出: [Lat, Lon, H] (deg, deg, m)
    Eigen::Vector3d GetBiasGyro() const;  // 输出: [Gx, Gy, Gz] (deg/h)
    Eigen::Vector3d GetBiasAcc() const;   // 输出: [Ax, Ay, Az] (ug)
    Eigen::Quaterniond GetQnb() const;    // 输出: 四元数

    void Step_Coarse(const Eigen::Vector3d& wm, const Eigen::Vector3d& vm);
    void Step_Fine(const Eigen::Vector3d& wm, const Eigen::Vector3d& vm);
    void Step_Nav(const Eigen::Vector3d& wm, const Eigen::Vector3d& vm);
    void Finish_Coarse();
    void Finish_Fine();    

    // ========== 混合对准接口 (推荐使用) ==========
    // 一键执行：粗对准 + KF精对准(姿态) + 几何均值(零偏)
    HybridAlignResult Run_HybridAlign(const std::vector<IMUData>& data, 
                                       const HybridAlignConfig& cfg = HybridAlignConfig());
    
    // 计算几何均值零偏 (给定姿态)
    Eigen::Vector3d ComputeGeometricEB(const Eigen::Vector3d& att,
                                        const std::vector<IMUData>& data,
                                        double eb_sigma_allan_dph = 0.003);
};

#endif // SINS_ENGINE_H