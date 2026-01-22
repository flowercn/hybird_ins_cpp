#pragma once
#include <Eigen/Dense>
#include <cmath>

using namespace Eigen;
 
class GLV {
public:
// 基础单位与转换常数
    double deg, rad, sec, min, hour, ppm;
    double dph, dpsh, arcmin, arcsec;

    // WGS-84 椭球体参数
    double Re, wie, f, e2, ge, g0;
    double m, beta, beta1;

    // 状态量 (遵循 pos = [lat, lon, h])
    Vector3d pos;
    double RMh, RNh;
    Vector3d wien;
    Vector3d gn;
    double g;

    // 传感器单位
    double ug, Hz, ugpsHz, ugpg2, nm;

    // 仿真与采样参数
    double tscale;
    double ts;
    int ns;

    // 补偿系数 (圆锥/划船)
    MatrixXd cs; 
    int csmax;
    int csCompensate;

    // 历史值 (上一时刻采样)
    Vector3d wm_1, vm_1;

    // 构造函数声明
    GLV();
};