#ifndef INSMATH_H
#define INSMATH_H

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "glv.h"

using namespace Eigen;
 
class INSMath {
public:
    // 1. 基础矢量运算
    static Matrix3d askew(const Vector3d& v);   // 反对称阵
    static Matrix3d rv2m(const Vector3d& rv);   // [新增] 旋转矢量转矩阵

    // 2. 姿态表示转换 (严格 3-1-2 顺序)
    static Matrix3d a2mat(const Vector3d& att);             // 欧拉角转矩阵
    static Quaterniond a2qua(const Vector3d& att);          // 欧拉角转四元数
    static Vector3d m2att(const Matrix3d& Cnb);             // 矩阵转欧拉角
    static Matrix3d q2mat(const Quaterniond& qnb);          // 四元数转矩阵
    static Quaterniond m2qua(const Matrix3d& Cnb);          // 矩阵转四元数
    static Vector3d q2rv(const Quaterniond& qnb);           // 四元数转旋转矢量

    // 3. 核心更新与评估算法
    static Quaterniond qupdt2(const Quaterniond& qnb0, const Vector3d& rv_ib, const Vector3d& rv_in); //
    static Vector3d aa2phi(const Vector3d& att1, const Vector3d& att0); // 计算失准角

    // 4. 辅助工具
    static Vector3d r2dms(double rad, const GLV& glv);      // 弧度转度分秒
    static Quaterniond rv2q(const Vector3d& rv);
};

#endif