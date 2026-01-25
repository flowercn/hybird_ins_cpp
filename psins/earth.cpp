#include "earth.h"
void Earth::eupdate(const Vector3d& pos, const Vector3d& vn) {
    // 纬度 pos(0), 经度 pos(1), 高度 pos(2)
    sl = sin(pos(0)); 
    cl = cos(pos(0)); 
    tl = sl / cl;
    sl2 = sl * sl;
    double rc = 1.0 - e2 * sl2;
    double sqrc = sqrt(rc);

    // 1. 地球曲率半径计算
    RMh = Re * (1.0 - e2) / (rc * sqrc) + pos(2);
    RNh = Re / sqrc + pos(2);
    clRNh = cl * RNh;

    // 2. 地球自转速率
    wien << 0.0, wie * cl, wie * sl;
    
    // 3. 牵连角速率 (vn: [vE, vN, vU])
    double vE_RNh = vn(0) / RNh;
    wenn << -vn(1) / RMh, vE_RNh, vE_RNh * tl;
    
    winn = wien + wenn;
    wcor = 2.0 * wien + wenn;

    // 4. 重力模型 (修正 s2L 变量名一致性)
    double s2L = 2.0 * sl * cl; 
    double gL = g0 * (1.0 + beta * sl2 - beta1 * (s2L * s2L));
    
    double hR = pos(2) / (Re * (1.0 - f * sl2));
    gL = gL * (1.0 - 2.0 * hR - 5.0 * hR * hR);
    gn << 0.0, 0.0, -gL;
    
    // 5. 补偿加速度 (Coriolis + Gravity)
    gcc = -wcor.cross(vn) + gn;

    // 存储中间变量
    this->sl = sl; this->cl = cl; this->tl = tl; this->sl2 = sl2;
}

// 在 Eth 类中添加构造函数
Earth::Earth(const GLV& glv) {
    // 基础常数同步
    this->Re = glv.Re; 
    this->e2 = glv.e2; 
    this->wie = glv.wie; 
    this->g0 = glv.g0;
    this->beta = glv.beta;
    this->beta1 = glv.beta1;
    this->f = glv.f;
    this->eupdate(glv.pos, Vector3d::Zero());
}