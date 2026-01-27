#include "sins_engine.h"
#include <iostream>

using namespace std;
using namespace Eigen;

SinsEngine::SinsEngine(double ts) 
    : eth(GLV()), 
      ins(Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero(), ts, eth) 
{
    ins.ts = ts;
}

void SinsEngine::Sins_Init(const Vector3d& init_pos, 
                           const Vector3d& init_vel, 
                           const Vector3d& init_att, 
                           const Vector3d& init_eb, 
                           const Vector3d& init_db) {
    // 1. 刷新地球模型
    eth = Earth(glv);
    eth.update(init_pos, init_vel);

    // 2. 填充初始结果
    res_init.valid = true;
    res_init.align_time = 0.0;
    
    res_init.pos = init_pos;
    res_init.vel = init_vel;
    res_init.att = init_att;
    res_init.eb  = init_eb;
    res_init.db  = init_db;

    // 3. 重置
    res_coarse.valid = false;
    res_fine.valid = false;
    
    coarse_sum_wm.setZero();
    coarse_sum_vm.setZero();
    coarse_timer = 0.0;
    
    cout << "[Sins] Initialized state directly." << endl;
}

void SinsEngine::SetKFConfig(const KFConfig& cfg) {
    this->kf_cfg = cfg;
}

// --- 粗对准 ---
void SinsEngine::Run_Coarse_Phase(const std::vector<IMUData>& data_chunk) {
    if (data_chunk.empty()) {
        cout << "[Coarse] Skipped (No Data)." << endl;
        return;
    }
    coarse_sum_wm.setZero();
    coarse_sum_vm.setZero();
    coarse_timer = 0.0;

    for (const auto& d : data_chunk) {
        Step_Coarse(d.wm, d.vm);
    }
    Finish_Coarse();
}

void SinsEngine::Step_Coarse(const Vector3d& wm, const Vector3d& vm) {
    coarse_sum_wm += wm;
    coarse_sum_vm += vm;
    coarse_timer += ins.ts;
}

void SinsEngine::Finish_Coarse() {
    if (coarse_timer < 0.1) return;

    Vector3d w_avg = coarse_sum_wm / coarse_timer;
    Vector3d f_avg = coarse_sum_vm / coarse_timer;
    
    // 解析对准
    // [Fix] 这里的 (0,0,1) 指的是当地地理坐标系的"天向"单位矢量，而非重力
    Vector3d up_ref(0, 0, 1); 
    double lat = res_init.pos(0);
    Vector3d wie_n(0, cos(lat), sin(lat));

    Vector3d vb = f_avg.normalized();
    Vector3d wb = w_avg.normalized();
    
    // 构造 n 系基准向量 (East, North, Up)
    Vector3d r_east_n = wie_n.cross(up_ref).normalized();
    Vector3d r_north_n = up_ref.cross(r_east_n).normalized(); // 注意叉乘顺序: Up x East = North
    Matrix3d Mn; Mn << r_east_n, r_north_n, up_ref;

    // 构造 b 系测量向量
    Vector3d r_east_b = wb.cross(vb).normalized();
    Vector3d r_north_b = vb.cross(r_east_b).normalized();
    Matrix3d Mb; Mb << r_east_b, r_north_b, vb;

    Matrix3d Cnb = Mn * Mb.transpose();
    
    res_coarse.valid = true;
    res_coarse.align_time = coarse_timer;
    res_coarse.att = INSMath::m2att(Cnb);
    res_coarse.vel.setZero();
    res_coarse.pos = res_init.pos;
    // 粗对准不估计零偏，清零以防万一
    res_coarse.eb.setZero();
    res_coarse.db.setZero();
}

// --- 精对准 (复用 INS 核心引擎) ---
void SinsEngine::Run_Fine_Phase(const std::vector<IMUData>& data_chunk) {
    if (data_chunk.empty()) {
        cout << "[Fine] Skipped (No Data)." << endl;
        return;
    }

    // 1. 确定初始状态 (继承粗对准或初始值)
    AlignResult* start = res_coarse.valid ? &res_coarse : &res_init;
    
    if (res_coarse.valid) cout << "[Fine] Starting with Coarse Alignment result." << endl;
    else cout << "[Fine] Starting with Init/Default result." << endl;

    // 2. 初始化 INS 引擎 (完全重置 INS 状态)
    // 这一步非常关键，把粗对准的姿态、初始速度(0)、位置、初始零偏(0或预设)装载进 INS
    eth.update(start->pos, start->vel);
    ins = INSState(start->att, start->vel, start->pos, ins.ts, eth);
    ins.set_bias(start->eb, start->db); // 设置初始零偏

    // 3. 初始化 KF
    // 注意：KF 的初始不确定度 (P阵) 依然来自 kf_cfg
    kf.Init(ins.ts, glv, 
            kf_cfg.phi_init_err(0), 
            kf_cfg.web_psd, 
            kf_cfg.wdb_psd, 
            kf_cfg.eb_sigma, 
            kf_cfg.db_sigma, 
            kf_cfg.wvn_err(0));

    // 4. 批量步进
    for (const auto& d : data_chunk) {
        Step_Fine(d.wm, d.vm);
    }

    Finish_Fine(); 
    res_fine.align_time = data_chunk.size() * ins.ts;
}

void SinsEngine::Step_Fine(const Vector3d& wm, const Vector3d& vm) {
    // 1. 调用高精度机械编排 (Don't Reinvent the Wheel!)
    // ins.update 内部会自动执行：
    //   a. 零偏补偿 (使用当前的 ins.eb/db)
    //   b. 圆锥/划船误差补偿
    //   c. 速度/位置更新 (含科氏力、重力)
    //   d. 姿态更新
    //   e. 地球参数刷新
    ins.update(wm, vm, glv, eth);

    // 2. 构造 KF 所需的输入量
    // KF 需要比力投影增量 dv_n (用于 F 阵中 f^n x phi 项)
    // 我们可以近似计算：Cnb * (vm - bias*ts)
    // 这里的 vm 最好是去过零偏的，虽然 KF 对 F 阵精度要求不高，但严谨点好
    Vector3d vm_pure = vm - ins.db * ins.ts; 
    Vector3d dv_n_measure = ins.Cnb * vm_pure; 

    // 3. KF 预测与更新
    // 观测值 Z = V_ins - V_ref (假设零速, V_ref=0) => Z = ins.vn
    kf.UpdatePhi(dv_n_measure, ins.Cnb, eth.wien, ins.ts);
    kf.Update(ins.vn); 
    
    // 4. [Fix] 全闭环反馈：直接修正 INS 对象的状态！
    // 传入引用：ins.qnb, ins.vn, ins.eb, ins.db, ins.Cnb
    // KF 内部会修改这些值，并清空 KF 自身的误差状态
    kf.Feedback(ins.qnb, ins.vn, ins.eb, ins.db, ins.Cnb);

    // 5. 状态同步与约束
    // 因为 Feedback 修改了 qnb，必须同步更新欧拉角 att，否则 GetAttDeg() 会拿旧数据
    ins.attsyn();
    
    // [工程经验] 静态对准期间，强制锁定位置和高度
    // 防止因加速度计零偏未收敛导致的位置发散，进而导致地球参数(g, wie)计算错误
    ins.pos = res_init.pos; 
    ins.vn(2) = 0.0; // 高程通道阻尼
}

void SinsEngine::Finish_Fine() {
    res_fine.valid = true;
    
    // 直接从 INS 引擎中提取最终状态
    // 因为是全闭环，ins.eb/db 已经是收敛后的最终零偏
    res_fine.att = ins.att;
    res_fine.vel = ins.vn;
    res_fine.pos = ins.pos;
    res_fine.eb = res_init.eb + kf.xk.segment<3>(6);
    res_fine.db = res_init.db + kf.xk.segment<3>(9);
    
    cout << "[Fine] Done. Bias Gyro: " << (res_fine.eb * glv.rad * 3600).transpose() << " deg/h" << endl;
    cout << "[Fine] Done. Bias Acc:  " << (res_fine.db / glv.ug).transpose() << " ug" << endl;
    Vector3d att_deg = res_fine.att * glv.rad;
    cout << "[Fine] Done. Att (deg):  P=" << att_deg(0) 
         << ", R=" << att_deg(1) 
         << ", Y=" << att_deg(2) << endl;
}

// --- 导航 (保持不变) ---
void SinsEngine::Run_Nav_Phase(const std::vector<IMUData>& data_chunk) {
    if (data_chunk.empty()) return;

    AlignResult* best = &res_init;
    if (res_fine.valid) {
        best = &res_fine;
        cout << "[Nav] Using FINE alignment result (" << best->align_time << "s)." << endl;
    }
    else if (res_coarse.valid) {
        best = &res_coarse;
        cout << "[Nav] Using COARSE alignment result (" << best->align_time << "s)." << endl;
    }
    else {
        cout << "[Nav] Using INITIAL default values (No alignment performed)." << endl;
    }

    eth.update(best->pos, best->vel);
    ins = INSState(best->att, best->vel, best->pos, ins.ts, eth);
    ins.set_bias(best->eb, best->db);

    for (const auto& d : data_chunk) {
        Step_Nav(d.wm, d.vm);
    }
    cout << "[Nav] Finished processing " << data_chunk.size() << " epochs." << endl;
}

void SinsEngine::Step_Nav(const Vector3d& wm, const Vector3d& vm) {
    ins.update(wm, vm, glv, eth);
    ins.vn(2) = 0.0; 
    ins.pos(2) = res_init.pos(2);
}

// Getters 
Vector3d SinsEngine::GetAttDeg() const { return ins.att * glv.rad; }
Vector3d SinsEngine::GetVel() const { return ins.vn; }
Vector3d SinsEngine::GetPosDeg() const {
    Vector3d p = ins.pos;
    p(0) *= glv.rad; p(1) *= glv.rad;
    return p;
}
Vector3d SinsEngine::GetBiasGyro() const { return ins.eb * glv.rad * 3600.0; }
Vector3d SinsEngine::GetBiasAcc() const { return ins.db / glv.ug; }
Quaterniond SinsEngine::GetQnb() const { return ins.qnb; }