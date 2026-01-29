#include "earth.h"
void Earth::update(const Vector3d& pos, const Vector3d& vn) {
    this->pos = pos;
    this->vn = vn;

    double lat = pos(0);
    double h = pos(2); 
    sl = sin(lat);
    cl = cos(lat);
    tl = sl / cl;
    sl2 = sl * sl;

    double rc = 1.0 - e2 * sl2;
    double sqrc = sqrt(rc);
    RMh = Re * (1.0 - e2) / (rc * sqrc) + h;
    RNh = Re / sqrc + h;
    clRNh = cl * RNh;

    wien << 0.0, wie * cl, wie * sl;
    double vE_RNh = vn(0) / RNh;
    wenn << -vn(1) / RMh, vE_RNh, vE_RNh * tl;
    winn = wien + wenn;

    double s2L = sin(2.0 * lat); 
    double gL = g0 * (1.0 + beta * sl2 - beta1 * (s2L * s2L));
    
    double hR = h / (Re * (1.0 - f * sl2));
    gL = gL * (1.0 - 2.0 * hR - 5.0 * hR * hR);

    gn << 0.0, 0.0, -gL;
    
    wcor = 2.0 * wien + wenn;
    gcc = -wcor.cross(vn) + gn;

    // 存储中间变量
}
