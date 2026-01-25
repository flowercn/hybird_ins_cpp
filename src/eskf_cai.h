#ifndef ESKF_CAI_H
#define ESKF_CAI_H

#include <Eigen/Dense>
#include <random>
#include <iostream>
#include "glv.h"
#include "ins_state.h"
#include "ins_math.h"
#include "earth.h" 

using namespace Eigen;

// ==========================================
// 1. 原子陀螺仿真器
// ==========================================
class CAIGSimulator {
public:
    double arw;       
    double bias_inst; 
    double T_cycle;   
    
    Vector3d bias_drift;     
    Vector3d w_ie_n_ref;     
    Matrix3d Cnb_anchor;     
    
    std::mt19937 gen;
    std::normal_distribution<double> norm_dist{0.0, 1.0};

    CAIGSimulator(double t_cycle_s = 2.0) : T_cycle(t_cycle_s) {
        double deg_to_rad = M_PI / 180.0;
        // MATLAB: 2e-4 deg/sqrt(h) -> 5e-4 (shot noise)
        // 这里使用你给出的原子参数
        double arw_input = 5.0e-4; 
        arw = arw_input * deg_to_rad / 60.0; 

        // bias instability
        double bi_input = 2.0e-5;
        bias_inst = bi_input * deg_to_rad / 3600.0;

        bias_drift.setZero();
    }

    void Init(const Quaterniond& qnb_aligned, double lat) {
        Cnb_anchor = INSMath::q2mat(qnb_aligned);
        GLV glv; 
        w_ie_n_ref << 0.0, glv.wie * cos(lat), glv.wie * sin(lat);
        std::random_device rd;
        gen.seed(rd());
    }

    Vector3d GetMeasurement() {
        Vector3d w_true_b = Cnb_anchor.transpose() * w_ie_n_ref;
        double sigma_white = arw / sqrt(T_cycle);
        Vector3d noise_white;
        noise_white << norm_dist(gen), norm_dist(gen), norm_dist(gen);
        noise_white *= sigma_white;

        // 简化的 bias instability 游走
        Vector3d noise_walk;
        noise_walk << norm_dist(gen), norm_dist(gen), norm_dist(gen);
        bias_drift += noise_walk * (bias_inst * sqrt(T_cycle)); 

        return w_true_b + bias_drift + noise_white;
    }
};

// ==========================================
// 2. 15维误差状态卡尔曼滤波 (MATLAB 复刻版)
// ==========================================
class ESKF15 {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    int n = 15;
    VectorXd xk;
    MatrixXd Pk;
    MatrixXd Qk;
    MatrixXd R_cai; 

    ESKF15() {
        xk = VectorXd::Zero(n);
        Pk = MatrixXd::Identity(n, n);
    }

    void Init(const INSState& ins, const GLV& glv) {
        xk.setZero();
        
        // Pk: 初始协方差 (对应 imuerr 中的初值或一般经验值)
        VectorXd P_diag(15);
        P_diag << 0.1*glv.arcmin, 0.1*glv.arcmin, 0.1*glv.arcmin, // phi (10 deg -> small)
                  0.1, 0.1, 0.1,                                  // dv
                  1.0, 1.0, 1.0,                                  // dpos
                  0.5*glv.dph, 0.5*glv.dph, 0.5*glv.dph,          // eb
                  100.0*glv.ug, 100.0*glv.ug, 100.0*glv.ug;       // db
        Pk = P_diag.array().square().matrix().asDiagonal();

        // Qk: 过程噪声 (对应 imuerr 中的 web/wdb)
        // MATLAB: web = 2e-4 deg/sqrt(h)
        // MATLAB: wdb = 5 ug/sqrt(Hz)
        
        double web_val = 2.0e-4 * glv.dpsh;   
        double wdb_val = 5.0 * glv.ugpsHz; 
        
        VectorXd Q_diag(15);
        Q_diag << VectorXd::Zero(9), 
                  Vector3d::Constant(web_val).array().square(), 
                  Vector3d::Constant(wdb_val).array().square();
        Qk = Q_diag.asDiagonal();
        Qk *= ins.ts;

        // Rk: CAI 观测噪声
        // MATLAB r_cai_val = 1.3225e-14
        R_cai = Matrix3d::Identity() * 1.3225e-14; 
    }

    // Predict: 严格复刻 MATLAB monitor 矩阵构建
    void Predict(const INSState& ins, const Earth& eth) {
        MatrixXd Ft = MatrixXd::Zero(n, n);
        
        Vector3d vn = ins.vn; 
        Matrix3d Cnb = ins.Cnb;
        Vector3d fn = ins.an; // 比力
        
        double f_RMh = 1.0 / eth.RMh;
        double f_RNh = 1.0 / eth.RNh;
        double f_clRNh = 1.0 / eth.clRNh;
        
        double vE_clRNh = vn(0) * f_clRNh;
        double vN_RMh2  = vn(1) * f_RMh * f_RMh; 
        double vE_RNh2  = vn(0) * f_RNh * f_RNh; 

        Matrix3d Map = Matrix3d::Zero();
        Map(1, 0) = -eth.wien(2); 
        Map(2, 0) = eth.wien(1) + vE_clRNh / eth.cl;
        Map(0, 2) = vN_RMh2;
        Map(1, 2) = -vE_RNh2;
        Map(2, 2) = -vE_RNh2 * eth.tl;

        Matrix3d Avn = INSMath::askew(vn);
        Matrix3d Maa = INSMath::askew(-eth.winn);

        Matrix3d Mav = Matrix3d::Zero();
        Mav(1, 0) = f_RNh;
        Mav(2, 0) = f_RNh * eth.tl;
        Mav(0, 1) = -f_RMh;

        Matrix3d Mva = INSMath::askew(-fn); 
        Matrix3d Mvv = Avn * Mav - INSMath::askew(eth.wcor);

        Matrix3d Map_cor = Map;
        Map_cor(1, 0) = Map(1, 0) - eth.wien(2); 
        Map_cor(2, 0) = Map(2, 0) + eth.wien(1); 
        
        Matrix3d Mvp = Avn * Map_cor;
        
        // 重力修正项
        double scl = eth.sl * eth.cl;
        double g0 = eth.g0; 
        Mvp(2, 0) = Mvp(2, 0) - g0 * (5.27094e-3 * 2 + 2.32718e-5 * 4 * eth.sl2) * scl;
        Mvp(2, 2) = Mvp(2, 2) + 3.086e-6;

        Matrix3d Mpv = ins.Mpv; 
        Matrix3d Mpp = Matrix3d::Zero();
        Mpp(1, 0) = vE_clRNh * eth.tl;
        Mpp(0, 2) = -vN_RMh2;
        Mpp(1, 2) = -vE_RNh2 / eth.cl;

        Ft.block<3,3>(0, 0) = Maa;
        Ft.block<3,3>(0, 3) = Mav;
        Ft.block<3,3>(0, 6) = Map;
        Ft.block<3,3>(0, 9) = -Cnb;

        Ft.block<3,3>(3, 0) = Mva;
        Ft.block<3,3>(3, 3) = Mvv;
        Ft.block<3,3>(3, 6) = Mvp;
        Ft.block<3,3>(3, 12) = Cnb;

        Ft.block<3,3>(6, 3) = Mpv;
        Ft.block<3,3>(6, 6) = Mpp;

        MatrixXd Phi = MatrixXd::Identity(n, n) + Ft * ins.ts;
        
        xk = Phi * xk;
        Pk = Phi * Pk * Phi.transpose() + Qk;
        Pk = 0.5 * (Pk + Pk.transpose());
    }

    // Update: 只保留 CAI (无 ZUPT)
    void Update(const Vector3d& w_cai_meas, const Vector3d& w_ins_avg) {
        Vector3d dw = w_cai_meas - w_ins_avg;
        
        VectorXd Z = dw; 
        MatrixXd H = MatrixXd::Zero(3, n);
        H.block<3,3>(0, 9) = -Matrix3d::Identity(); 
        
        MatrixXd R = R_cai; 
        
        MatrixXd PHt = Pk * H.transpose();
        MatrixXd S = H * PHt + R;
        MatrixXd K = PHt * S.inverse();
        
        xk = xk + K * (Z - H * xk);
        MatrixXd I = MatrixXd::Identity(n, n);
        Pk = (I - K * H) * Pk * (I - K * H).transpose() + K * R * K.transpose();
    }

    void Feedback(INSState& ins) {
        Vector3d phi = xk.segment<3>(0);
        Quaterniond q_corr = INSMath::rv2q(phi);
        ins.qnb = q_corr * ins.qnb;
        ins.attsyn(); 
        
        ins.vn  -= xk.segment<3>(3);
        ins.pos -= xk.segment<3>(6);
        
        ins.eb += xk.segment<3>(9);
        ins.db += xk.segment<3>(12);
        
        xk.setZero();
    }
};

#endif