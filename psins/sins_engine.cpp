#include "sins_engine.h"
#include "ins_math.h"
#include <iomanip>
#include <cmath>

using namespace std;
using namespace Eigen;

SinsEngine::SinsEngine(double ts) : state(EngineState::IDLE), coarse_duration(0), fine_duration(0) {
    ins.ts = ts;
    align_finished = false;
    kf_round = 0;
    scale_ratio = 1.0;
    accum_bias_gyro.setZero();
    accum_bias_acc.setZero();
}

void SinsEngine::Init(const AlignConfig& cfg, double coarse_time_s, double fine_time_s) {
    this->align_cfg = cfg;
    this->coarse_duration = coarse_time_s;
    this->fine_duration = fine_time_s;
    
    // 设定第二轮 KF 的切换时间点 (例如精对准进行到一半时)
    // 这里设定为: 粗对准结束后的 100 秒进行切换
    this->time_kf_switch = coarse_duration + 100.0; 
    if (this->time_kf_switch >= coarse_duration + fine_duration) {
        this->time_kf_switch = coarse_duration + fine_duration * 0.5;
    }

    eth = Earth(glv);
    
    coarse_acc_wm.setZero();
    coarse_acc_vm.setZero();
    coarse_sample_cnt = 0;
    kf_round = 0;
    scale_ratio = 1.0;
    accum_bias_gyro.setZero();
    accum_bias_acc.setZero();

    state = EngineState::COARSE_ALIGNING;
    
    cout << "------------------------------------------------" << endl;
    cout << "[SinsEngine] Initialized Robust Iterative Alignment:" << endl;
    cout << "  Phase 1: Analytic Coarse Align (0 ~ " << coarse_duration << " s)" << endl;
    cout << "  Phase 2: KF Round 1 (Large P0) (" << coarse_duration << " ~ " << time_kf_switch << " s)" << endl;
    cout << "  Phase 3: KF Round 2 (Refinement) (" << time_kf_switch << " ~ " << (coarse_duration + fine_duration) << " s)" << endl;
    cout << "------------------------------------------------" << endl;
}

// 静态粗对准逻辑 (保持不变)
Vector3d SinsEngine::AlignCoarse(const Vector3d& w_avg, const Vector3d& f_avg, double latitude) {
    Vector3d gn_ref(0, 0, 1); 
    Vector3d wie_n(0, std::cos(latitude), std::sin(latitude));
    Vector3d fb_meas = f_avg.normalized();
    Vector3d wb_meas = w_avg.normalized();

    // 构造导航系矩阵 Mn
    Vector3d r_up_n = gn_ref; 
    Vector3d r_east_n = wie_n.cross(r_up_n).normalized(); 
    Vector3d r_north_n = r_up_n.cross(r_east_n).normalized();
    Matrix3d Mn; Mn << r_east_n, r_north_n, r_up_n;

    // 构造机体系矩阵 Mb
    Vector3d r_up_b = fb_meas; 
    Vector3d r_east_b = wb_meas.cross(r_up_b).normalized();
    Vector3d r_north_b = r_up_b.cross(r_east_b).normalized();
    Matrix3d Mb; Mb << r_east_b, r_north_b, r_up_b;

    Matrix3d Cnb = Mn * Mb.transpose();
    return INSMath::m2att(Cnb);
}

void SinsEngine::Step(const Vector3d& wm, const Vector3d& vm, double t) {
    if (wm.hasNaN() || vm.hasNaN()) {
        cerr << "[SINS-FATAL] NaN input detected at t=" << t << endl;
        state = EngineState::FAULT;
        return;
    }
    if (state == EngineState::FAULT) return;

    // ==============================================================================
    // 阶段 1: 粗对准 (Accumulate & Average)
    // ==============================================================================
    if (state == EngineState::COARSE_ALIGNING) {
        coarse_acc_wm += wm; 
        coarse_acc_vm += vm; 
        coarse_sample_cnt++;

        if (t >= coarse_duration) {
            cout << "\n[SinsEngine] >>> Coarse Alignment Finished at t=" << t << "s" << endl;

            // 1. 计算均值
            double total_time = coarse_sample_cnt * ins.ts;
            Vector3d w_avg = coarse_acc_wm / total_time; 
            Vector3d f_avg = coarse_acc_vm / total_time; 

            // 2. 自动计算 Scale Ratio (解决重力模长不匹配问题)
            eth.eupdate(align_cfg.pos_ref, Vector3d::Zero());
            double local_g = eth.gn.norm();
            double data_g  = f_avg.norm();
            this->scale_ratio = local_g / data_g;
            
            cout << "  [Auto-Scale] Local G: " << local_g << ", Data G: " << data_g 
                 << ", Ratio: " << scale_ratio << endl;
            
            // 应用 Scale
            f_avg *= scale_ratio;

            // 3. 执行解析粗对准
            // 使用算出来的 att_coarse，而不是配置的 att_ref！
            Vector3d att_coarse = SinsEngine::AlignCoarse(w_avg, f_avg, align_cfg.pos_ref(0));
            cout << "  [Result] Coarse Att: " << (att_coarse * glv.rad).transpose() << " (Will be used as Initial Guess)" << endl;

            // 4. 初始化 KF Round 1
            // 策略：使用 Coarse 结果初始化，但给非常大的 P0，允许 KF 大幅修正
            eth.eupdate(align_cfg.pos_ref, Vector3d::Zero());
            align_vn.setZero();
            align_pos = align_cfg.pos_ref;
            align_qnb = INSMath::a2qua(att_coarse); // <--- 关键：使用粗对准结果
            
            align_nn = 2; 
            align_nts = align_nn * ins.ts;
            align_step_counter = 0;
            
            // 初始化 KF
            align_kf.Init(align_nts, glv, 0.0, align_cfg.web_psd, align_cfg.wdb_psd, align_cfg.eb_sigma, align_cfg.db_sigma, align_cfg.wvn_err(0));    

            // 设置宽松的 P0 (Large Uncertainty)
            // 例如：水平 1度，方位 5度
            Vector3d large_phi_err; 
            large_phi_err << 1.0, 1.0, 5.0; 
            large_phi_err *= glv.deg;
            align_kf.Pxk.block<3,3>(0,0) = large_phi_err.array().square().matrix().asDiagonal();
            
            kf_round = 1;
            state = EngineState::FINE_ALIGNING;
            cout << "[SinsEngine] Switched to KF Round 1 (Pull-in Phase, Large P0)." << endl;
        }
    }
    
    // ==============================================================================
    // 阶段 2: 精对准 (Iterative KF)
    // ==============================================================================
    else if (state == EngineState::FINE_ALIGNING) {
        
        // ------------------------------------------------------------------
        // 数据预处理：应用 Scale 和 上一轮固化的 Bias
        // ------------------------------------------------------------------
        Vector3d vm_scaled = vm * this->scale_ratio;
        Vector3d wm_comp = wm - accum_bias_gyro * ins.ts;
        Vector3d vm_comp = vm_scaled - accum_bias_acc * ins.ts;

        // ------------------------------------------------------------------
        // 双子样累加
        // ------------------------------------------------------------------
        if (align_step_counter == 0) {
            align_last_data.wm = wm_comp;
            align_last_data.vm = vm_comp;
            align_last_data.t = t;
            align_step_counter++;
            return; // 等待下一帧
        }
        
        // align_step_counter == 1
        const auto& d1 = align_last_data;
        Vector3d wmm = d1.wm + wm_comp;
        Vector3d vmm = d1.vm + vm_comp;
        
        // 误差补偿
        Vector3d dphim = 2.0/3.0 * d1.wm.cross(wm_comp);
        Vector3d phim = wmm + dphim; 
        Vector3d scullm = 2.0/3.0 * (d1.wm.cross(vm_comp) + d1.vm.cross(wm_comp));
        Vector3d rotm = 0.5 * wmm.cross(vmm);
        Vector3d dvbm = vmm + rotm + scullm;

        // 机械编排
        Matrix3d Cnn = INSMath::rv2m(-eth.wien * align_nts / 2.0); 
        Matrix3d Cnb = INSMath::q2mat(align_qnb); 
        Vector3d dvn = Cnn * Cnb * dvbm;
        align_vn = align_vn + dvn + eth.gn * align_nts;
        
        Vector3d rv_phi = phim;
        Vector3d rv_zeta = eth.winn * align_nts; 
        align_qnb = INSMath::qupdt2(align_qnb, rv_phi, rv_zeta);

        // KF Update
        align_kf.UpdatePhi(dvn, Cnb, eth.wien, align_nts);
        align_kf.Update(align_vn); 
        align_kf.Feedback(align_qnb, align_vn);
        
        align_step_counter = 0;

        // ------------------------------------------------------------------
        // 逻辑检查：切换 Round 2 或 结束
        // ------------------------------------------------------------------
        
        // [切换到 Round 2]
        if (kf_round == 1 && t >= time_kf_switch) {
            cout << "\n[SinsEngine] >>> KF Round 1 Finished at t=" << t << "s." << endl;
            
            // 1. 提取当前 KF 估计出的零偏
            Vector3d est_eb = align_kf.xk.segment<3>(6);
            Vector3d est_db = align_kf.xk.segment<3>(9);
            
            // 2. 固化到 accum_bias 中
            accum_bias_gyro += est_eb;
            accum_bias_acc  += est_db;
            
            cout << "  [Round 1 Result] EB: " << (accum_bias_gyro * glv.rad * 3600.0).transpose() << " deg/h" << endl;
            cout << "  [Round 1 Result] DB: " << (accum_bias_acc / glv.ug).transpose() << " ug" << endl;
            
            // 3. 重置 KF (Re-Init)
            align_kf.xk.setZero(); 
            align_kf.Pxk.setIdentity(); 
            
            // 4. 设置精细的 P0 (Small Uncertainty)
            Vector3d small_phi_err; 
            small_phi_err << 0.02, 0.02, 0.1; // 0.1度航向误差
            small_phi_err *= glv.deg;
            
            // ========================= [修正点] =========================
            // 原错误代码: Vector3d p_diag(12); 
            // 修正为 VectorXd (动态大小) 或 Matrix<double, 12, 1>
            VectorXd p_diag(12); 
            // ============================================================

            p_diag << small_phi_err(0), small_phi_err(1), small_phi_err(2), // phi
                      0.1, 0.1, 0.1,                                        // dv
                      align_cfg.eb_sigma * 0.5, align_cfg.eb_sigma * 0.5, align_cfg.eb_sigma * 0.5, // eb
                      align_cfg.db_sigma * 0.5, align_cfg.db_sigma * 0.5, align_cfg.db_sigma * 0.5; // db
            
            align_kf.Pxk = p_diag.array().square().matrix().asDiagonal();
            
            kf_round = 2;
            cout << "[SinsEngine] Switched to KF Round 2 (Refinement, Small P0, Bias Fixed)." << endl;
        }
        
        // [精对准完全结束]
        if (t >= (coarse_duration + fine_duration)) {
             cout << "\n[SinsEngine] >>> Fine Alignment (All Rounds) Finished." << endl;
             
             // 最终结果 = 系统状态 + KF残差状态
             AlignResult res;
             res.att = INSMath::m2att(INSMath::q2mat(align_qnb)); 
             res.pos = align_pos;
             res.vn  = align_vn;
             
             // 最终零偏 = 固化零偏 + KF残差零偏
             res.eb = accum_bias_gyro + align_kf.xk.segment<3>(6);
             res.db = accum_bias_acc  + align_kf.xk.segment<3>(9);
             
             cout << "  [Nav Start] Att:   " << (res.att * glv.rad).transpose() << endl;
             cout << "  [Nav Start] EB:    " << (res.eb * glv.rad * 3600.0).transpose() << endl;
             cout << "  [Nav Start] DB:    " << (res.db / glv.ug).transpose() << endl;
             
             // 初始化导航状态
             eth.eupdate(align_cfg.pos_ref, Vector3d::Zero());
             ins = INSState(res.att, Vector3d::Zero(), align_cfg.pos_ref, ins.ts, eth);
             ins.set_bias(res.eb, res.db); 
             
             ins.is_align = false;
             ins.vertical_damping_mode = 1; 
             ins.vn(2) = 0.0;               

             state = EngineState::NAVIGATING; 
        }
    }
    
    // ==============================================================================
    // 阶段 3: 导航
    // ==============================================================================
    else if (state == EngineState::NAVIGATING) {
        // 注意：导航阶段也要应用 scale_ratio ！
        Vector3d vm_scaled = vm * this->scale_ratio;
        ins.update(wm, vm_scaled, glv, eth);
        
        if (ins.qnb.coeffs().hasNaN() || ins.pos.hasNaN()) {
            cerr << "[SINS-FATAL] Algorithm Diverged at t=" << t << endl;
            state = EngineState::FAULT;
        }
    }
}

// Getters ... (保持不变)
Vector3d SinsEngine::GetAttDeg() const {
    if (state == EngineState::FINE_ALIGNING) {
        return INSMath::m2att(INSMath::q2mat(align_qnb)) * glv.rad;
    } 
    else if (state == EngineState::COARSE_ALIGNING) {
        return Vector3d::Zero();
    }
    return ins.att * glv.rad;
}

Vector3d SinsEngine::GetVel() const {
    if (state == EngineState::FINE_ALIGNING) return align_vn;
    return ins.vn; 
}

Vector3d SinsEngine::GetPosDeg() const {
    Vector3d p;
    if (state == EngineState::FINE_ALIGNING) p = align_pos;
    else if (state == EngineState::NAVIGATING) p = ins.pos;
    else p = align_cfg.pos_ref; 

    p(0) *= glv.rad; p(1) *= glv.rad;
    return p;
}

Vector3d SinsEngine::GetBiasGyro() const {
    // 返回总零偏
    if (state == EngineState::FINE_ALIGNING) return accum_bias_gyro + align_kf.xk.segment<3>(6);
    if (state == EngineState::NAVIGATING) return ins.eb;
    return Vector3d::Zero();
}

Vector3d SinsEngine::GetBiasAcc() const {
    if (state == EngineState::FINE_ALIGNING) return accum_bias_acc + align_kf.xk.segment<3>(9);
    if (state == EngineState::NAVIGATING) return ins.db;
    return Vector3d::Zero();
}

void SinsEngine::InjectBias(const Vector3d& db_gyro, const Vector3d& db_acc) {
    if (state == EngineState::NAVIGATING) {
        ins.set_bias(ins.eb + db_gyro, ins.db + db_acc);
    }
}

Quaterniond SinsEngine::GetQnb() const {
    if (state == EngineState::FINE_ALIGNING) return align_qnb;
    if (state == EngineState::NAVIGATING) return ins.qnb;
    return Quaterniond::Identity(); 
}