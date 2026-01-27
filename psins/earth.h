#pragma once
#include <Eigen/Dense>
#include "glv.h"

using namespace Eigen;

class Earth {
public:
    double Re, wie, f, e2, g0;
    double beta, beta1;

    Vector3d pos;  // 当前位置 (lat, lon, h)
    Vector3d vn;   // 当前速度 (vn, ve, vd)
    
    double sl, cl, tl, sl2; // sin(lat), cos(lat), tan(lat)...
    double RMh, RNh, clRNh; // 子午圈/卯酉圈半径
    
    Vector3d wien; // 地球自转角速度 (n系)
    Vector3d wenn; // 运输率 (n系)
    Vector3d winn; // wien + wenn
    Vector3d gn;   // 当地重力向量
    Vector3d gcc;  // 哥氏/向心加速度 (有害加速度)

    Earth() = default;
    Earth(const GLV& glv) {
        this->Re = glv.Re;
        this->wie = glv.wie;
        this->f = glv.f;
        this->e2 = glv.e2;
        this->g0 = glv.g0;
        this->beta = glv.beta;
        this->beta1 = glv.beta1;
    }

    void update(const Vector3d& pos, const Vector3d& vn);
};