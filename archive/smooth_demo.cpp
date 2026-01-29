/**
 * @file smooth_demo.cpp
 * @brief 回溯平滑算法演示 - 原子陀螺辅助纯惯导
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "sins_engine.h"
#include "support.h"
#include "cai_sim.h"

using namespace std;
using namespace Eigen;

namespace OptimalParams {
    const Vector3d att_deg(0.0344, 0.335, 0.590);
    const Vector3d eb_dph(-0.0017, -0.0026, -0.0103);
}

int main() {
    cout << "==================================================" << endl;
    cout << "  Backtracking Smoothing - Correct Implementation" << endl;
    cout << "==================================================" << endl;
    
    string csv_path = "../fog3h.csv";
    double ts = 1.0 / 400.0;
    GLV glv;
    
    Vector3d pos_nj(32.0286 * glv.deg, 118.8533 * glv.deg, 17.0);
    
    Earth eth(glv);
    eth.update(pos_nj, Vector3d::Zero());
    double local_g = eth.gn.norm();
    
    cout << "\n[Optimal Parameters (Ground Truth)]" << endl;
    cout << "  Attitude (deg): " << OptimalParams::att_deg.transpose() << endl;
    cout << "  Bias (deg/h):   " << OptimalParams::eb_dph.transpose() << endl;
    
    Vector3d att_true = OptimalParams::att_deg * glv.deg;
    Vector3d eb_true = OptimalParams::eb_dph * glv.deg / 3600.0;
    
    cout << "\n[1] Loading data..." << endl;
    auto imu_data = LoadIMUData(csv_path, ts, local_g);
    cout << "    Samples: " << imu_data.size() << " (" << imu_data.size() * ts / 60 << " min)" << endl;
    
    cout << "\n[2] FOG Scale Factor Calibration..." << endl;
    
    Vector3d wm_sum = Vector3d::Zero();
    for (const auto& d : imu_data) wm_sum += d.wm;
    Vector3d w_fog_avg = wm_sum / (imu_data.size() * ts);
    double fog_norm = w_fog_avg.norm();
    double scale_factor = glv.wie / fog_norm;
    
    cout << "    FOG measured |w| = " << fixed << setprecision(4) << fog_norm / glv.deg * 3600 << " deg/h" << endl;
    cout << "    Scale factor K   = " << setprecision(6) << scale_factor << endl;
    
    for (auto& d : imu_data) d.wm *= scale_factor;
    
    wm_sum.setZero();
    for (const auto& d : imu_data) wm_sum += d.wm;
    w_fog_avg = wm_sum / (imu_data.size() * ts);
    cout << "    FOG avg (deg/h): " << (w_fog_avg / glv.deg * 3600).transpose() << endl;
    
    SinsEngine engine(ts);
    engine.Sins_Init(pos_nj, Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero());
    
    // 使用统一的原子陀螺仿真器
    cai::AtomicGyroSimulator atom(pos_nj, glv);
    atom.Init(att_true);
    atom.PrintConfig();
    
    double nav_time = 3600.0;
    size_t nav_samples = min(static_cast<size_t>(nav_time / ts), imu_data.size());
    const double T_CYCLE = atom.GetParams().T_cycle;
    const int SAMPLES_PER_CYCLE = static_cast<int>(T_CYCLE / ts);
    
    // ========== 场景1: eb = eb_true ==========
    cout << "\n[3] Scenario 1: eb = eb_true..." << endl;
    engine.eth.update(pos_nj, Vector3d::Zero());
    engine.ins = INSState(att_true, Vector3d::Zero(), pos_nj, ts, engine.eth);
    engine.ins.set_bias(eb_true, Vector3d::Zero());
    
    double max_drift_1 = 0;
    for (size_t i = 0; i < nav_samples; ++i) {
        engine.Step_Nav(imu_data[i].wm, imu_data[i].vm);
        double dN = (engine.ins.pos(0) - pos_nj(0)) * glv.Re;
        double dE = (engine.ins.pos(1) - pos_nj(1)) * cos(pos_nj(0)) * glv.Re;
        double drift = sqrt(dN*dN + dE*dE);
        if (drift > max_drift_1) max_drift_1 = drift;
    }
    cout << "    Drift: " << fixed << setprecision(3) << max_drift_1 / 1852.0 << " nm/h" << endl;
    
    // ========== 场景2: eb = 0 ==========
    cout << "\n[4] Scenario 2: eb = 0 (no compensation)..." << endl;
    engine.eth.update(pos_nj, Vector3d::Zero());
    engine.ins = INSState(att_true, Vector3d::Zero(), pos_nj, ts, engine.eth);
    engine.ins.set_bias(Vector3d::Zero(), Vector3d::Zero());
    
    double max_drift_2 = 0;
    for (size_t i = 0; i < nav_samples; ++i) {
        engine.Step_Nav(imu_data[i].wm, imu_data[i].vm);
        double dN = (engine.ins.pos(0) - pos_nj(0)) * glv.Re;
        double dE = (engine.ins.pos(1) - pos_nj(1)) * cos(pos_nj(0)) * glv.Re;
        double drift = sqrt(dN*dN + dE*dE);
        if (drift > max_drift_2) max_drift_2 = drift;
    }
    cout << "    Drift: " << fixed << setprecision(3) << max_drift_2 / 1852.0 << " nm/h" << endl;
    
    // ========== 场景3: 回溯校正 ==========
    cout << "\n[5] Scenario 3: Backtracking with atom gyro..." << endl;
    
    // 正确的 wie_b 参考 (原子陀螺测量)
    Vector3d wie_b_ref = atom.GetTrueWieB();
    cout << "    wie_b_ref (deg/h): " << (wie_b_ref / glv.deg * 3600).transpose() << endl;
    
    // 第一遍：收集每周期的 eb
    vector<Vector3d> eb_per_cycle;
    Vector3d fog_sum = Vector3d::Zero();
    
    for (size_t i = 0; i < nav_samples; ++i) {
        fog_sum += imu_data[i].wm;
        if ((i + 1) % SAMPLES_PER_CYCLE == 0) {
            Vector3d eb_this = (fog_sum - wie_b_ref * T_CYCLE) / T_CYCLE;
            eb_per_cycle.push_back(eb_this);
            fog_sum.setZero();
        }
    }
    
    // 第二遍：回放
    engine.eth.update(pos_nj, Vector3d::Zero());
    engine.ins = INSState(att_true, Vector3d::Zero(), pos_nj, ts, engine.eth);
    
    double max_drift_3 = 0;
    size_t cycle_idx = 0;
    
    for (size_t i = 0; i < nav_samples; ++i) {
        if (cycle_idx < eb_per_cycle.size()) {
            engine.ins.eb = eb_per_cycle[cycle_idx];
        }
        engine.Step_Nav(imu_data[i].wm, imu_data[i].vm);
        if ((i + 1) % SAMPLES_PER_CYCLE == 0) cycle_idx++;
        
        double dN = (engine.ins.pos(0) - pos_nj(0)) * glv.Re;
        double dE = (engine.ins.pos(1) - pos_nj(1)) * cos(pos_nj(0)) * glv.Re;
        double drift = sqrt(dN*dN + dE*dE);
        if (drift > max_drift_3) max_drift_3 = drift;
    }
    cout << "    Drift: " << fixed << setprecision(3) << max_drift_3 / 1852.0 << " nm/h" << endl;
    
    Vector3d eb_avg = Vector3d::Zero();
    for (const auto& eb : eb_per_cycle) eb_avg += eb;
    eb_avg /= eb_per_cycle.size();
    cout << "    eb avg (deg/h): " << (eb_avg / glv.deg * 3600).transpose() << endl;
    
    // ==========================================================
    cout << "\n==================================================" << endl;
    cout << "  Results Summary" << endl;
    cout << "==================================================" << endl;
    cout << fixed << setprecision(3);
    cout << "  [1] eb = eb_true:     " << max_drift_1 / 1852.0 << " nm/h" << endl;
    cout << "  [2] eb = 0:           " << max_drift_2 / 1852.0 << " nm/h" << endl;
    cout << "  [3] Backtracking:     " << max_drift_3 / 1852.0 << " nm/h" << endl;
    cout << "==================================================" << endl;
    
    return 0;
}
