#include "ins_state.h"
#include "ins_math.h"
#include <iostream>

INSState::INSState(const Vector3d& att0, const Vector3d& vn0, const Vector3d& pos0, double sampling_ts, const Earth& eth) {
    this->ts = sampling_ts;
    this->att = att0; 
    this->vn = vn0; 
    this->pos = pos0;
    
    // 备份初始状态
    this->vn0 = vn0; 
    this->pos0 = pos0;

    // 初始化姿态矩阵
    this->Cnb = INSMath::a2mat(att0);
    this->Cnb0 = this->Cnb;
    this->qnb = INSMath::a2qua(att0); 
    
    // 初始化零偏
    this->eb.setZero(); 
    this->db.setZero();
    this->an.setZero();
    this->Kg.setZero();
    this->Ka.setZero();

    // 初始化 Mpv (位置更新矩阵)
    this->Mpv.setZero();
    this->Mpv(0, 1) = 1.0 / eth.RMh;
    this->Mpv(1, 0) = 1.0 / eth.clRNh;
    this->Mpv(2, 2) = 1.0;

    // 【关键】初始化历史数据为 0 (接管 GLV 的工作)
    this->wm_last.setZero(); 
    this->vm_last.setZero();
    
}

void INSState::set_bias(const Vector3d& new_eb, const Vector3d& new_db) {
    this->eb = new_eb;
    this->db = new_db;
}

void INSState::set_scalefactor(const Vector3d& new_kg, const Vector3d& new_ka) {
    this->Kg = new_kg;
    this->Ka = new_ka;
}

// 圆锥划船补偿
void INSState::cnscl(const Vector3d& wm, const Vector3d& vm, Vector3d& phim, Vector3d& dvbm) {
    // dphim = 1/12 * (w_{k-1} x w_k)
    Vector3d dphim = (1.0/12.0) * wm_last.cross(wm);
    // scullm = 1/12 * (w_{k-1} x v_k + v_{k-1} x w_k)
    Vector3d scullm = (1.0/12.0) * (wm_last.cross(vm) + vm_last.cross(wm));
    Vector3d rotm = 0.5 * wm.cross(vm);
    
    phim = wm + dphim;
    dvbm = vm + rotm + scullm;
    
    wm_last = wm; 
    vm_last = vm;
}

void INSState::update(const Vector3d& wm, const Vector3d& vm, const GLV& glv, Earth& eth) {
    double nts2 = ts / 2.0;

    //this->debug_vm_raw = vm;

    //static int debug_cnt = 0;
    // 1. 扣除零偏 (Bias Correction)
    Vector3d wm_unbiased = wm - eb * ts;
    Vector3d vm_unbiased = vm - db * ts;

    /// 2. 扣除刻度因数误差 (Scale Factor Correction)
    // 推荐写法：转为 array() 进行逐元素运算，然后再转回 matrix()
    Vector3d wm_pure = wm_unbiased.array() / (Vector3d::Ones() + Kg).array();
    Vector3d vm_pure = vm_unbiased.array() / (Vector3d::Ones() + Ka).array();

    //this->debug_vm_pure = vm_pure;
    // 2. 补偿 (计算圆锥/划船效应)
    Vector3d phim, dvbm;
    cnscl(wm_pure, vm_pure, phim, dvbm);

    // 3. 速度/位置外推 (用于计算中间时刻的地球参数)
    Vector3d vn_mid = vn + an * nts2; 
    Vector3d pos_mid = pos + (Mpv * vn_mid) * nts2;
    
    // 更新 Earth 参数 (使用外推后的位置)
    eth.update(pos_mid, vn_mid); 

    // 4. 速度更新
    Vector3d dv_n = Cnb * dvbm; 
    Vector3d vn1 = vn + dv_n + eth.gcc * ts;
    an = dv_n / ts + eth.gcc; // 记录加速度用于下一次外推

    // 5. 位置更新 (梯形积分)
    Mpv(0, 1) = 1.0 / eth.RMh;
    Mpv(1, 0) = 1.0 / eth.clRNh;
    Vector3d pos1 = pos + (Mpv * (vn + vn1) * 0.5) * ts;

    // 6. 姿态更新
    // 使用 winn (n系转动) 修正姿态
    qnb = INSMath::qupdt2(qnb, phim, eth.winn * ts);
    attsyn(); 

    // 7. 更新状态
    vn = vn1;
    pos = pos1;
}

void INSState::attsyn() {
    qnb.normalize();
    Cnb = INSMath::q2mat(qnb);
    att = INSMath::m2att(Cnb);
}