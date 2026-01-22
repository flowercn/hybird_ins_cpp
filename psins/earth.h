#pragma once
#include <Eigen/Dense>
#include "glv.h" // 引用物理常数类

using namespace Eigen;

class Earth {
public:
    // 物理常数备份 (由 GLV 初始化时提供)
    double Re, wie, f, e2, g0, beta, beta1;

    // 实时计算的地理参数 (n系)
    double RMh, RNh, clRNh;
    Vector3d wien, wenn, winn, wcor;
    Vector3d gn, gcc;

    // 中间变量存储
    double sl, cl, tl, sl2;

    // 默认构造
    Earth() = default;
    // 对应 MATLAB 的 einit(): 使用 GLV 常数初始化
    Earth(const GLV& glv);
    void eupdate(const Vector3d& pos, const Vector3d& vn);
}; 