#include "ins_state.h"
#include "ins_math.h"
#include <iostream>

INSState::INSState() 
    : ts(0.01), is_open_loop(false), is_align(false), csCompensate(1), vertical_damping_mode(0) {
    att.setZero(); vn.setZero(); pos.setZero(); 
    an.setZero(); eb.setZero(); db.setZero(); 
    Cnb.setIdentity(); qnb.setIdentity(); Mpv.setZero();
    wm_1.setZero(); vm_1.setZero(); phim.setZero(); dvbm.setZero();
    vn0.setZero(); pos0.setZero(); Cnb0.setIdentity();
}

INSState::INSState(const Vector3d& att0, const Vector3d& vn0, const Vector3d& pos0, double sampling_ts, const Earth& eth) {
    this->ts = sampling_ts;
    this->att = att0; 
    this->vn = vn0; 
    this->pos = pos0;
    
    this->vn0 = vn0; 
    this->pos0 = pos0;

    this->Cnb = INSMath::a2mat(att0);
    this->qnb = INSMath::a2qua(att0); 
    this->Cnb0 = this->Cnb;

    this->an.setZero();
    this->eb.setZero(); 
    this->db.setZero();
    
    this->Mpv.setZero();
    this->Mpv(0, 1) = 1.0 / eth.RMh;
    this->Mpv(1, 0) = 1.0 / eth.clRNh;
    this->Mpv(2, 2) = 1.0;

    this->wm_1.setZero(); 
    this->vm_1.setZero();
    
    this->is_open_loop = false; 
    this->is_align = false; 
    this->csCompensate = 1;
    this->vertical_damping_mode = 1; 
}

INSState::INSState(const VectorXd& avp0, double sampling_ts, const Earth& eth) 
    : INSState(avp0.segment<3>(0), avp0.segment<3>(3), avp0.segment<3>(6), sampling_ts, eth) {}

void INSState::set_bias(const Vector3d& new_eb, const Vector3d& new_db) {
    this->eb = new_eb;
    this->db = new_db;
}

void INSState::cnscl(const Vector3d& wm, const Vector3d& vm) {
    if (csCompensate == 0) {
        this->phim = wm;
        this->dvbm = vm;
        return;
    }
    Vector3d dphim = (1.0/12.0) * wm_1.cross(wm);
    Vector3d scullm = (1.0/12.0) * (wm_1.cross(vm) + vm_1.cross(wm));
    Vector3d rotm = 0.5 * wm.cross(vm);
    
    this->phim = wm + dphim;
    this->dvbm = vm + rotm + scullm;
    
    wm_1 = wm; 
    vm_1 = vm;
}

void INSState::update(const Vector3d& wm, const Vector3d& vm, const GLV& glv, Earth& eth) {
    double nts2 = ts / 2.0;

    // [Step 1] 零偏补偿
    Vector3d wm_pure = wm - eb * ts;
    Vector3d vm_pure = vm - db * ts;

    // [Step 2] 圆锥/划船
    this->cnscl(wm_pure, vm_pure);

    // [Step 3] 外推
    Vector3d vn_mid = vn + an * nts2; 
    Vector3d pos_mid = pos + (Mpv * vn_mid) * nts2;
    eth.eupdate(pos_mid, vn_mid); 

    // [Step 4] 速度更新
    Vector3d dv_n = Cnb * dvbm; 
    Vector3d vn1 = vn + dv_n + eth.gcc * ts;
    an = dv_n / ts + eth.gcc; 

    // [Step 5] 位置更新
    Mpv(0, 1) = 1.0 / eth.RMh;
    Mpv(1, 0) = 1.0 / eth.clRNh;
    Vector3d pos1 = pos + (Mpv * (vn + vn1) * 0.5) * ts;

    // [Step 6] 姿态更新
    qnb = INSMath::qupdt2(qnb, this->phim, eth.winn * ts);
    attsyn(); 

    // [Step 7] 迭代
    vn = vn1;
    pos = pos1;

    // [Step 8] 约束与阻尼
    if (is_open_loop) {
        vn.setZero();
        pos = pos0;
    } else if (vertical_damping_mode == 1) {
        vn(2) = 0.0;       
        pos(2) = pos0(2);  
    } else if (is_align) {
    }
}

void INSState::attsyn() {
    qnb.normalize();
    Cnb = INSMath::q2mat(qnb);
    att = INSMath::m2att(Cnb);
}