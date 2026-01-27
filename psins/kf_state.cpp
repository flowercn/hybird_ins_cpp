#include "kf_state.h"
#include "ins_math.h"
#include <iostream>

KFAlignVN::KFAlignVN() {
    xk.setZero(n);
    Pxk.setIdentity(n, n);
    Phikk_1.setIdentity(n, n);
    Hk.setZero(m, n);
    Hk.block<3,3>(0, 3) = Matrix3d::Identity();
}

void KFAlignVN::Init(double nts, const GLV& glv, 
                     double phi0, double web_psd, double wdb_psd, 
                     double eb_sigma, double db_sigma, double wvn_err) {
    
    // --- Qk (过程噪声) ---
    Qk.setZero(n, n);
    Vector3d web2 = Vector3d::Constant(web_psd).array().square(); 
    Vector3d wdb2 = Vector3d::Constant(wdb_psd).array().square(); 
    
    Qk.diagonal().segment<3>(0) = web2 * nts;
    Qk.diagonal().segment<3>(3) = wdb2 * nts; // 速度随机游走
    Qk.diagonal().segment<3>(6) = Vector3d::Zero(); // 零偏通常建模为随机常数(0噪声)或一阶马尔可夫
    Qk.diagonal().segment<3>(9) = Vector3d::Zero(); 
    
    // --- Rk (测量噪声) ---
    Rk = Matrix3d::Identity() * (wvn_err * wvn_err) / nts;
    
    // --- Pxk (初始协方差) ---
    VectorXd p_diag(12);
    p_diag << phi0, phi0, phi0,                 // phi
              0.1, 0.1, 0.1,                    // dv
              eb_sigma, eb_sigma, eb_sigma,     // eb
              db_sigma, db_sigma, db_sigma;     // db
    Pxk = p_diag.array().square().matrix().asDiagonal();
    Phikk_1.setIdentity();
}

void KFAlignVN::UpdatePhi(const Vector3d& dvn, const Matrix3d& Cnb, const Vector3d& wnie, double nts) {
    Phikk_1.setIdentity(); 
    
    // 1. 姿态部分: -w_in x
    Matrix3d Ft_att = INSMath::askew(-wnie);
    Phikk_1.block<3,3>(0,0) += Ft_att * nts;
    
    // 2. 速度误差项: 
    // d(dv)/dt = f x phi - (2w_ie + w_en) x dv
    // f x phi 部分:
    Phikk_1.block<3,3>(3,0) = INSMath::askew(dvn) * nts; // 注意：必须乘 nts (离散化)
    
    // [Fix] 补充科氏力项 -(2w_ie + w_en) x dv
    // 这里的 wnie = w_ie + w_en。近似认为 2*w_ie + w_en ≈ 2*wnie (对于低速对准误差很小)
    // 否则需要单独传入 wie。此处用 askew(-2.0 * wnie) 近似
    Phikk_1.block<3,3>(3,3) += INSMath::askew(-2.0 * wnie) * nts;
    
    // 3. 零偏项: 
    Matrix3d Cnbts = Cnb * nts;
    // phi 对 eb: -Cnb * dt (注意符号)
    Phikk_1.block<3,3>(0,6) = -Cnbts;
    
    // v 对 db: +Cnb * dt
    Phikk_1.block<3,3>(3,9) = Cnbts;
}

void KFAlignVN::Update(const Vector3d& vn_meas) {
    // 预测
    xk = Phikk_1 * xk;
    Pxk = Phikk_1 * Pxk * Phikk_1.transpose() + Qk;

    // 更新
    // 注意: vn_meas 必须是 (v_ins - v_ref) 或者是 (v_ins) 且 Hx=v_ins?
    // 标准写法: Z = v_ins - v_ref. H=[0 I 0 0]. Z_pred = Hx = dv.
    // Innovation yk = Z - Z_pred = (v_ins - v_ref) - dv.
    VectorXd Hx = Hk * xk;
    yk = vn_meas - Hx; 
    
    MatrixXd S = Hk * Pxk * Hk.transpose() + Rk;
    Kk = Pxk * Hk.transpose() * S.inverse();
    
    xk = xk + Kk * yk;
    MatrixXd I = MatrixXd::Identity(n, n);
    Pxk = (I - Kk * Hk) * Pxk;
}

void KFAlignVN::Feedback(Quaterniond& qnb, Vector3d& vn, Vector3d& eb, Vector3d& db, Matrix3d& Cnb) {
    double fb_ratio = 0.90; // 反馈比例
    double remain = 1.0 - fb_ratio;

    // 1. 姿态反馈
    Vector3d phi_fb = xk.segment<3>(0) * fb_ratio;
    // 修正 qnb (根据 dv_dot = f x phi 推导，修正为正旋转)
    Quaterniond q_corr = INSMath::rv2q(phi_fb); 
    qnb = q_corr * qnb;  
    qnb.normalize();
    
    // [Fix] 立即同步 Cnb，防止 INS 下一步计算比力投影时出错
    Cnb = INSMath::q2mat(qnb);

    xk.segment<3>(0) *= remain;
    
    // 2. 速度反馈
    Vector3d vn_fb = xk.segment<3>(3) * fb_ratio;
    vn -= vn_fb;
    xk.segment<3>(3) *= remain;

    // 3. [Fix] 零偏反馈
    // 陀螺零偏 (eb): 对应 Phikk_1 中的 -Cnb，暗示 x_eb = (b_hat - b_true)
    // 所以我们需要减去估计出的误差
    Vector3d eb_fb = xk.segment<3>(6) * fb_ratio;
    eb -= eb_fb; 
    xk.segment<3>(6) *= remain;

    // 加计零偏 (db): 对应 Phikk_1 中的 +Cnb，暗示 x_db = (b_true - b_hat)
    // 所以我们需要加上估计出的误差
    Vector3d db_fb = xk.segment<3>(9) * fb_ratio;
    db += db_fb;
    xk.segment<3>(9) *= remain;
}