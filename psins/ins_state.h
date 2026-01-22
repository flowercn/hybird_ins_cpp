#ifndef INS_STATE_H
#define INS_STATE_H

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "glv.h"   
#include "earth.h" 

using namespace Eigen;
 
class INSState {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // --- 状态量 ---
    Vector3d att;      // 欧拉角 (rad)
    Vector3d vn;       // 速度 (m/s)
    Vector3d pos;      // 位置 (rad, rad, m)
    Vector3d an;       // 导航系加速度 (m/s^2)
    Vector3d eb, db;   // 零偏 (rad/s, m/s^2)
    
    Matrix3d Cnb;      
    Matrix3d Mpv;      
    Quaterniond qnb;   

    // --- 备份与配置 ---
    Vector3d vn0, pos0;
    Matrix3d Cnb0;
    double ts;             
    
    // --- 模式开关 ---
    bool is_open_loop;     // 开环模式 (vn=0)
    bool is_align;         // 对准模式 (通常无阻尼)
    int csCompensate;      // 1:开启圆锥划船补偿
    int vertical_damping_mode; 

    // --- 中间变量 ---
    Vector3d phim, dvbm;   
    Vector3d wm_1, vm_1;   

    // --- 构造函数 ---
    INSState(); 
    INSState(const Vector3d& att0, const Vector3d& vn0, const Vector3d& pos0, double sampling_ts, const Earth& eth);
    INSState(const VectorXd& avp0, double sampling_ts, const Earth& eth);

    // --- 核心方法 ---
    void update(const Vector3d& wm, const Vector3d& vm, const GLV& glv, Earth& eth);
    void cnscl(const Vector3d& wm, const Vector3d& vm);
    void attsyn();
    void set_bias(const Vector3d& new_eb, const Vector3d& new_db);
};

#endif // INS_STATE_H