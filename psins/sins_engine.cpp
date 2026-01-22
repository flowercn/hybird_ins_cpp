#include "sins_engine.h"
#include <iomanip>

using namespace std;
using namespace Eigen;

SinsEngine::SinsEngine(double ts) : state(EngineState::IDLE), align_duration(0.0) {
    ins.ts = ts;
}

void SinsEngine::Init(const AlignConfig& cfg, double align_time_sec) {
    this->align_cfg = cfg;
    this->align_duration = align_time_sec;
    eth = Earth(glv);
    aligner.Init(ins.ts, glv, cfg);
    state = EngineState::ALIGNING;
    cout << "[SinsEngine] Initialized. Target Alignment Time: " << align_duration << " s." << endl;
}

void SinsEngine::Step(const Vector3d& wm, const Vector3d& vm, double t) {
    if (wm.hasNaN() || vm.hasNaN()) {
        cerr << "[SINS-FATAL] NaN input detected at t=" << t << endl;
        state = EngineState::FAULT;
        return;
    }
    if (state == EngineState::FAULT) return;
    IMUData curr_data;
    curr_data.wm = wm;
    curr_data.vm = vm;
    curr_data.t = t;
    if (state == EngineState::ALIGNING) {
        aligner.Step(curr_data);
        if (t >= align_duration) {
            cout << "\n[SinsEngine] >>> Time Reached (" << t << "s). Switching to Navigation..." << endl;
            AlignResult res = aligner.GetResult();
            if (!res.success) {
                cerr << "[SINS-ERROR] Alignment not converged/failed!" << endl;
                state = EngineState::FAULT;
                return;
            }

            // 2. 初始化纯惯导状态 (INSState)
            // 注意：这里使用对准算出的 att，和配置中的参考 pos (通常对准不更新位置或位置由GNSS给定)
            eth.eupdate(align_cfg.pos_ref, Vector3d::Zero());
            
            ins = INSState(res.att, Vector3d::Zero(), align_cfg.pos_ref, ins.ts, eth);
            ins.set_bias(res.eb, res.db); // 继承估计出的零偏
            ins.is_align = false;
            ins.vertical_damping_mode = 1; // 开启高度阻尼
            ins.vn(2) = 0.0;
            // 3. 切换状态
            state = EngineState::NAVIGATING;
            
            // 打印诊断信息
            cout << "  [Transition] Initial Attitude (deg): " << (res.att * glv.rad).transpose() << endl;
            cout << "  [Transition] Est Gyro Bias (deg/h):  " << (res.eb * glv.rad * 3600.0).transpose() << endl;
            cout << "  [Transition] Est Acc Bias (ug):      " << (res.db / glv.ug).transpose() << endl;
        }
    }
    // Case B: 正在导航阶段
    else if (state == EngineState::NAVIGATING) {
        
        // 执行机械编排
        ins.update(wm, vm, glv, eth);

        // [可选] 如果需要严格锁定高度，可以解除下方注释：
        // ins.pos(2) = align_cfg.pos_ref(2);

        // 输出保护 (Output Check)
        if (ins.qnb.coeffs().hasNaN() || ins.pos.hasNaN()) {
            cerr << "[SINS-FATAL] Algorithm Diverged at t=" << t << endl;
            state = EngineState::FAULT;
        }
    }
}

// --- Getter 实现 ---

Vector3d SinsEngine::GetAttDeg() const {
    if (state == EngineState::ALIGNING) {
        // 对准阶段返回 aligner 内部姿态
        return INSMath::m2att(INSMath::q2mat(aligner.qnb)) * glv.rad;
    }
    return ins.att * glv.rad;
}

Vector3d SinsEngine::GetVel() const {
    if (state == EngineState::ALIGNING) return aligner.vn;
    return ins.vn;
}

Vector3d SinsEngine::GetPosDeg() const {
    Vector3d p = (state == EngineState::ALIGNING) ? aligner.pos : ins.pos;
    p(0) *= glv.rad; p(1) *= glv.rad;
    return p;
}

Vector3d SinsEngine::GetBiasGyro() const {
    return (state == EngineState::ALIGNING) ? aligner.kf.xk.segment<3>(6) : ins.eb;
}

Vector3d SinsEngine::GetBiasAcc() const {
    return (state == EngineState::ALIGNING) ? aligner.kf.xk.segment<3>(9) : ins.db;
}

void SinsEngine::InjectBias(const Vector3d& db_gyro, const Vector3d& db_acc) {
    if (state == EngineState::NAVIGATING) {
        ins.set_bias(ins.eb + db_gyro, ins.db + db_acc);
    }
}

Quaterniond SinsEngine::GetQnb() const {
    if (state == EngineState::ALIGNING) {
        return aligner.qnb; // 需要确认 AlignmentEngine 定义了 qnb 为 public
    }
    return ins.qnb;
}