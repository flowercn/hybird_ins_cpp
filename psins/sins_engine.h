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
    Eigen::Vector3d kg;          // 陀螺刻度误差
    Eigen::Vector3d ka;          // 加计刻度误差
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

// 混合对准结果
struct HybridAlignResult {
    bool valid = false;
    Eigen::Vector3d att;     // 最终姿态 (rad)
    Eigen::Vector3d eb;      // 最终几何均值陀螺零偏 (rad/s)
    Eigen::Vector3d db;      // 最终几何均值加计零偏 (m/s^2) [NEW: 现在这里会存入有效值]
    double align_time;       // 对准总时长
    
    // --- 中间结果 (用于调试/日志打印) ---
    Eigen::Vector3d att_coarse;  // 粗对准姿态
    Eigen::Vector3d att_kf;      // KF 精对准结束时的姿态
    
    // 陀螺中间变量
    Eigen::Vector3d eb_raw;      // 原始测量残差 (rad/s)
    Eigen::Vector3d eb_scale;    // 陀螺几何约束缩放因子 (0.0 ~ 1.0)
    
    // [NEW] 加计中间变量 (对应代码中的 db_raw 和 scale_db)
    Eigen::Vector3d db_raw;      // 原始加计测量残差 (m/s^2)
    Eigen::Vector3d db_scale;    // 加计几何约束缩放因子 (0.0 ~ 1.0)
};

// 混合对准配置
struct HybridAlignConfig {
    double t_coarse = 60.0;         // 粗对准时长 (s)
    double t_fine = 3600.0;         // 精对准时长 (s) 
    
    // 零偏稳定性约束 (用于几何均值滤波)
    double eb_sigma_allan = 0.003;  // 陀螺 Allan 零偏稳定性 (deg/h)
    
    // [NEW] 加计稳定性约束
    // 建议设为 50.0 ~ 100.0 ug (基于您的 Allan 方差图底噪为 1ug，这就足够抓住 57ug 的误差了)
    double db_sigma_allan = 50.0;   // 加计 Allan 零偏稳定性 (ug)
    
    bool verbose = true;            // 是否输出过程信息
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
                   const Eigen::Vector3d& init_db,
                   const Eigen::Vector3d& init_kg, 
                   const Eigen::Vector3d& init_ka);
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