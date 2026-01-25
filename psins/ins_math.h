#ifndef INSMATH_H
#define INSMATH_H

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "glv.h"

using namespace Eigen;

class INSMath {
public:
    // 基础运算
    static Matrix3d askew(const Vector3d& v);   // 反对称矩阵 [v×]

    // 姿态转换 (严格 3-1-2 顺序)
    static Matrix3d rv2m(const Vector3d& rv);               // RV -> Cnb
    static Quaterniond rv2q(const Vector3d& rv);            // RV -> qnb
    static Matrix3d a2mat(const Vector3d& att);             // Euler -> Cnb
    static Quaterniond a2qua(const Vector3d& att);          // Euler -> qnb
    static Vector3d m2att(const Matrix3d& Cnb);             // Cnb -> Euler
    static Matrix3d q2mat(const Quaterniond& qnb);          // qnb -> Cnb
    static Quaterniond m2qua(const Matrix3d& Cnb);          // Cnb -> qnb
    static Vector3d q2rv(const Quaterniond& qnb);           // qnb -> RV

    // 核心算法
    static Quaterniond qupdt2(const Quaterniond& qnb0, const Vector3d& rv_ib, const Vector3d& rv_in); 
    static Vector3d aa2phi(const Vector3d& att1, const Vector3d& att0); // 计算失准角 phi

    // 辅助工具
    static Vector3d r2dms(double rad, const GLV& glv);      // Rad -> DMS
};

#endif