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
    Vector3d eb;       // 陀螺零偏 (rad/s)
    Vector3d db;       // 加计零偏 (m/s^2)
    Vector3d Kg;       // 陀螺刻度误差
    Vector3d Ka;       // 加计刻度误差
    
    Matrix3d Cnb;      
    Matrix3d Mpv;      
    Quaterniond qnb;   

    Vector3d wm_last;
    Vector3d vm_last;

    //Vector3d debug_vm_raw;  // 原始速度增量
    //Vector3d debug_vm_pure; // 修正后的速度增量

    Vector3d vn0, pos0;
    Matrix3d Cnb0;
    double ts;             
    
    // --- 构造函数 ---
    INSState() = delete;
    INSState(const Vector3d& att0, const Vector3d& vn0, const Vector3d& pos0, double ts, const Earth& eth);
    // --- core ---
    void set_bias(const Vector3d& new_eb, const Vector3d& new_db);
    void set_scalefactor(const Vector3d& new_kg, const Vector3d& new_ka);
    void update(const Vector3d& wm, const Vector3d& vm, const GLV& glv, Earth& eth);
    void attsyn();
private:
    void cnscl(const Vector3d& wm, const Vector3d& vm, Vector3d& phim, Vector3d& dvbm);
};

#endif // INS_STATE_H