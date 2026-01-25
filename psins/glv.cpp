#include "glv.h"

GLV::GLV() {
    // 基础单位与常数
    deg = M_PI / 180.0; //deg->rad
    rad = 180.0 / M_PI; //rad->deg
    sec = 1.0;
    min = 60.0;
    hour = 3600.0;
    ppm = 1e-6;
    dph = deg / hour; // 1dph=pi/180/3600 rad/s
    dpsh = deg / sqrt(hour); // 1dpsh=pi/180/sqrt(3600) rad/sqrt(s)
    arcmin = deg / 60.0; // 1arcmin=pi/180/60 rad
    arcsec = deg / 3600.0; // 1arcsec=pi/180/3600 rad

    // WGS-84 椭球体参数
    Re = 6.378136998405e6;
    wie = 7.2921151467e-5;
    f = 1.0 / 298.257223563;
    e2 = 2.0 * f - f * f;
    ge = 9.780325333434361;
    g0 = 9.7803267715;

    // 重力参数计算
    m = (wie * wie * Re) / ge;
    beta = 2.5 * m - f;
    beta1 = 0.125 * (5.0 * m * f - f * f);

    // 初始位置设置：pos = [lat, lon, h] (rad, rad, m)
    pos << 32.028611111111111 * deg, 118.85333333333333e+02 * deg, 17.0;

    // 初始地球参数计算 (einit 逻辑)
    double sl = sin(pos(0)), cl = cos(pos(0));
    double sl2 = sl * sl;
    double rc = 1.0 - e2 * sl2;
    double sqrc = sqrt(rc);

    RMh = Re * (1.0 - e2) / (rc * sqrc) + pos(2);
    RNh = Re / sqrc + pos(2);
    
    // 地球自转向量 (n系)
    wien << 0.0, wie * cl, wie * sl;

    // 初始重力计算
    double s2L = 2.0 * sl * cl;
    double gL = g0 * (1.0 + beta * sl2 - beta1 * (s2L * s2L));
    double hR = pos(2) / (Re * (1.0 - f * sl2));
    gL = gL * (1.0 - 2.0 * hR - 5.0 * hR * hR);
    
    gn << 0.0, 0.0, -gL;
    g = std::abs(gn(2)); 
    
    // 其他单位
    ug = 1e-6 * g;
    Hz = 1.0;
    ugpsHz = ug / sqrt(Hz);
    ugpg2 = ug / (g * g);
    nm = Re * arcmin;

    // 仿真/采样设置
    tscale = 1.0;
    ns = 1;
    ts = 1.0 / 400.0;

    // 补偿系数初始化 (圆锥/划船补偿)
    cs.setZero(5, 5); 
    cs.row(0) << 2.0/3.0,    0.0,      0.0,      0.0,      0.0;
    cs.row(1) << 9.0/20.0,   27.0/20.0, 0.0,      0.0,      0.0;
    cs.row(2) << 54.0/105.0, 92.0/105.0, 214.0/105.0, 0.0,    0.0;
    cs.row(3) << 250.0/504.0, 525.0/504.0, 650.0/504.0, 1375.0/504.0, 0.0;
    cs.row(4) << 2315.0/4620.0, 4558.0/4620.0, 7296.0/4620.0, 7834.0/4620.0, 15797.0/4620.0;
    
    csmax = 6; // size(glv.cs, 1) + 1
    csCompensate = 0;
}