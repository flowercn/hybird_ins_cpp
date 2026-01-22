#include "kf_state.h"
#include "ins_math.h"
#include <iostream>

KFAlignVN::KFAlignVN() {
    xk.setZero(n);
    Pxk.setIdentity(n, n);
    Phikk_1.setIdentity(n, n);
    Hk.setZero(m, n);
    // kf.Hk = [zeros(3),eye(3),zeros(3,6)]; 
    // MATLAB indices 4:6 -> C++ indices 3:5
    Hk.block<3,3>(0, 3) = Matrix3d::Identity();
}

void KFAlignVN::Init(double nts, const GLV& glv, 
                     double phi0, double web_psd, double wdb_psd, 
                     double eb_sigma, double db_sigma, double wvn_err) {
    
    // --- Qk (过程噪声) ---
    // MATLAB: kf.Qk = diag([imuerr.web; imuerr.wdb; zeros(6,1)])^2*nts;
    Qk.setZero(n, n);
    Vector3d web2 = Vector3d::Constant(web_psd).array().square(); // (rad/s)^2 ? or psd
    Vector3d wdb2 = Vector3d::Constant(wdb_psd).array().square(); 
    
    Qk.diagonal().segment<3>(0) = web2 * nts;
    Qk.diagonal().segment<3>(3) = wdb2 * nts;
    
    // --- Rk (测量噪声) ---
    // MATLAB: kf.Rk = diag(wvn)^2/nts;  <-- 注意这里除以 nts
    Rk = Matrix3d::Identity() * (wvn_err * wvn_err) / nts;
    
    // --- Pxk (初始协方差) ---
    // MATLAB: kf.Pxk = diag([phi0; [1;1;1]; imuerr.eb; imuerr.db])^2;
    VectorXd p_diag(12);
    p_diag << phi0, phi0, phi0,                 // phi
              1.0, 1.0, 1.0,                    // dv
              eb_sigma, eb_sigma, eb_sigma,     // eb
              db_sigma, db_sigma, db_sigma;     // db
    Pxk = p_diag.array().square().matrix().asDiagonal();
    
    // 初始 Phikk_1 固定部分: eye(12) + Ft*nts
    // Ft(1:3,1:3) = askew(-wnie)
    Phikk_1.setIdentity();
    // 动态部分在 UpdatePhi 中处理，这里只初始化
}

void KFAlignVN::UpdatePhi(const Vector3d& dvn, const Matrix3d& Cnb, const Vector3d& wnie, double nts) {
    // 复刻 alignvn.m 循环内的 Phikk_1 更新
    // kf.Phikk_1 = eye(12) + Ft*nts 基础
    // 动态更新部分：
    
    Phikk_1.setIdentity(); 
    
    // 1. 姿态部分: askew(-wnie)*nts
    Matrix3d Ft_att = INSMath::askew(-wnie);
    Phikk_1.block<3,3>(0,0) += Ft_att * nts;
    
    // 2. 速度误差项: kf.Phikk_1(4:6,1:3) = askew(dvn);
    Phikk_1.block<3,3>(3,0) = INSMath::askew(dvn);
    
    // 3. 零偏项: 
    // kf.Phikk_1(1:3,7:9) = -Cnbts;
    Matrix3d Cnbts = Cnb * nts;
    Phikk_1.block<3,3>(0,6) = -Cnbts;
    
    // kf.Phikk_1(4:6,10:12) = Cnbts;
    Phikk_1.block<3,3>(3,9) = Cnbts;
}

void KFAlignVN::Update(const Vector3d& vn_meas) {
    // 标准 KF 流程
    xk = Phikk_1 * xk;
    Pxk = Phikk_1 * Pxk * Phikk_1.transpose() + Qk;

    VectorXd Hx = Hk * xk;
    yk = vn_meas - Hx; 
    
    MatrixXd S = Hk * Pxk * Hk.transpose() + Rk;
    Kk = Pxk * Hk.transpose() * S.inverse();
    
    xk = xk + Kk * yk;
    MatrixXd I = MatrixXd::Identity(n, n);
    Pxk = (I - Kk * Hk) * Pxk;
}

void KFAlignVN::Feedback(Quaterniond& qnb, Vector3d& vn) {
    // 1. 姿态反馈
    // 对应 MATLAB qdelphi: qnb = qmul(rv2q(phi), qpb) -> 左乘正 phi
    Vector3d phi_fb = xk.segment<3>(0) * 0.91;
    
    // [回滚修正] 使用正 phi_fb，因为我们现在补全了 rotm，物理模型与 MATLAB 一致了
    Quaterniond q_err = INSMath::rv2q(phi_fb); 
    
    qnb = q_err * qnb;  
    qnb.normalize();
    xk.segment<3>(0) *= 0.09;
    
    // 2. 速度反馈 (保持不变)
    Vector3d vn_fb = xk.segment<3>(3) * 0.91;
    vn -= vn_fb;
    xk.segment<3>(3) *= 0.09;
}