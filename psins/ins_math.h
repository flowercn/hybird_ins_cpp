#ifndef INSMATH_H
#define INSMATH_H

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "glv.h"

using namespace Eigen;

namespace INSMath {

    // --- 基础运算 ---
    Matrix3d askew(const Vector3d& v);

    // --- 旋转矢量 <-> 矩阵/四元数 ---
    Matrix3d rv2m(const Vector3d& rv);
    Quaterniond rv2q(const Vector3d& rv);

    // --- 欧拉角 <-> 矩阵/四元数 (用户自定义顺序) ---
    Matrix3d a2mat(const Vector3d& att);
    Quaterniond a2qua(const Vector3d& att);
    
    // --- 四元数 <-> 矩阵/旋转矢量 ---
    Matrix3d q2mat(const Quaterniond& qnb);
    Vector3d q2rv(const Quaterniond& qnb);

    // --- 矩阵 <-> 四元数/欧拉角 ---
    Quaterniond m2qua(const Matrix3d& Cnb);
    Vector3d m2att(const Matrix3d& Cnb);
    
    // --- 核心更新 ---
    Quaterniond qupdt2(const Quaterniond& qnb0, const Vector3d& rv_ib, const Vector3d& rv_in);
    Vector3d aa2phi(const Vector3d& att1, const Vector3d& att0);

    // --- 工具 ---
    Vector3d r2dms(double rad, const GLV& glv);
};

#endif