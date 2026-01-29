/**
 * @file strategy_compare.cpp
 * @brief 多种回溯平滑策略对比分析
 * 
 * 关键发现：
 * - FOG Z轴角速率在1小时内从 7.959 下降到 ~7.8 deg/h（持续漂移）
 * - 原子陀螺保持稳定在 ~7.97 deg/h
 * - 这说明FOG存在低频零偏漂移，原子陀螺可以提供稳定基准
 * 
 * 策略：
 * 1. None: 不修正（基线）
 * 2. First Window: 用第一个窗口的 eb_diff
 * 3. Per Window: 每个窗口用各自的 eb_diff（跟踪漂移）
 * 4. Global Mean: 用所有窗口的 eb_diff 均值
 * 5. Linear Fit: 线性拟合 eb(t)
 * 6. Sliding Average: 滑动平均
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <random>
#include <cmath>
#include "sins_engine.h"
#include "support.h"

using namespace std;
using namespace Eigen;

// =============================================================================
// 最优参数
// =============================================================================
namespace OptimalParams {
    const Vector3d att_deg(0.0344, 0.335, 0.590);
    const Vector3d eb_dph(-0.0017, -0.0026, -0.0103);
}

// =============================================================================
// 原子陀螺仿真器
// =============================================================================
class AtomicGyroSimulator {
public:
    static constexpr double T_CYCLE = 2.0;
    static constexpr double ARW_DPSH = 2.0e-4;
    static constexpr double BIAS_DPH = 1.0e-5;
    
    double angle_noise_rad;
    Vector3d bias_vec;
    Vector3d wie_b;
    
    mt19937 rng;
    normal_distribution<double> noise_dist;
    
    AtomicGyroSimulator(const Vector3d& att_true_rad, const Vector3d& pos_rad, const GLV& glv)
        : rng(42)
    {
        double lat = pos_rad(0);
        Vector3d wie_n(0, glv.wie * cos(lat), glv.wie * sin(lat));
        Matrix3d Cnb = INSMath::a2mat(att_true_rad);
        wie_b = Cnb.transpose() * wie_n;
        
        double t_cycle_hours = T_CYCLE / 3600.0;
        double angle_noise_deg = ARW_DPSH * sqrt(t_cycle_hours);
        angle_noise_rad = angle_noise_deg * glv.deg;
        
        noise_dist = normal_distribution<double>(0.0, 1.0);
        
        double bias_rad_s = (BIAS_DPH * glv.deg) / 3600.0;
        bias_vec = Vector3d(bias_rad_s, bias_rad_s, bias_rad_s);
    }
    
    Vector3d Measure() { return wie_b + bias_vec; }
    Vector3d GetAngleNoise() {
        return Vector3d(noise_dist(rng), noise_dist(rng), noise_dist(rng)) * angle_noise_rad;
    }
    Vector3d GetTrueWieB() const { return wie_b; }
};

// =============================================================================
// 窗口数据
// =============================================================================
struct WindowData {
    double start_time;
    double end_time;
    Vector3d w_fog_mean;    // FOG均值 (rad/s)
    Vector3d w_atom_mean;   // 原子均值 (rad/s)
    Vector3d eb_diff;       // FOG-Atom ≈ eb_fog (rad/s)
};

// =============================================================================
// 线性拟合
// =============================================================================
void LinearFit(const vector<double>& t, const vector<double>& y, double& a, double& b) {
    int n = t.size();
    double sum_t = 0, sum_y = 0, sum_tt = 0, sum_ty = 0;
    for (int i = 0; i < n; ++i) {
        sum_t += t[i]; sum_y += y[i];
        sum_tt += t[i] * t[i]; sum_ty += t[i] * y[i];
    }
    double denom = n * sum_tt - sum_t * sum_t;
    if (fabs(denom) < 1e-12) { a = sum_y / n; b = 0; return; }
    b = (n * sum_ty - sum_t * sum_y) / denom;
    a = (sum_y - b * sum_t) / n;
}

// =============================================================================
// 收集窗口数据
// =============================================================================
vector<WindowData> CollectWindowData(
    const vector<IMUData>& imu_data,
    AtomicGyroSimulator& atom,
    double ts,
    double nav_time_s,
    double T_WINDOW,
    const GLV& glv)
{
    size_t nav_samples = min(static_cast<size_t>(nav_time_s / ts), imu_data.size());
    const int ATOM_SAMPLES_PER_CYCLE = static_cast<int>(AtomicGyroSimulator::T_CYCLE / ts);
    const int SAMPLES_PER_WINDOW = static_cast<int>(T_WINDOW / ts);
    
    atom.rng.seed(42);
    
    vector<WindowData> windows;
    Vector3d fog_sum = Vector3d::Zero();
    Vector3d atom_sum = Vector3d::Zero();
    double window_start = 0;
    int sample_in_window = 0;
    
    for (size_t i = 0; i < nav_samples; ++i) {
        fog_sum += imu_data[i].wm;
        
        if ((i + 1) % ATOM_SAMPLES_PER_CYCLE == 0) {
            Vector3d w_atom = atom.Measure();
            Vector3d dtheta_atom = w_atom * AtomicGyroSimulator::T_CYCLE + atom.GetAngleNoise();
            atom_sum += dtheta_atom;
        }
        
        sample_in_window++;
        
        if (sample_in_window >= SAMPLES_PER_WINDOW) {
            WindowData wd;
            wd.start_time = window_start;
            wd.end_time = (i + 1) * ts;
            wd.w_fog_mean = fog_sum / T_WINDOW;
            wd.w_atom_mean = atom_sum / T_WINDOW;
            wd.eb_diff = wd.w_fog_mean - wd.w_atom_mean;
            windows.push_back(wd);
            
            window_start = wd.end_time;
            fog_sum.setZero();
            atom_sum.setZero();
            sample_in_window = 0;
        }
    }
    
    return windows;
}

// =============================================================================
// 使用指定零偏进行导航
// =============================================================================
double NavigateWithBias(
    SinsEngine& engine,
    const vector<IMUData>& imu_data,
    const Vector3d& att_init,
    const Vector3d& pos_ref,
    const vector<Vector3d>& eb_per_window,
    double ts,
    double T_WINDOW,
    double nav_time_s,
    const GLV& glv)
{
    size_t nav_samples = min(static_cast<size_t>(nav_time_s / ts), imu_data.size());
    const int SAMPLES_PER_WINDOW = static_cast<int>(T_WINDOW / ts);
    
    engine.eth.update(pos_ref, Vector3d::Zero());
    engine.ins = INSState(att_init, Vector3d::Zero(), pos_ref, ts, engine.eth);
    
    double max_drift = 0;
    size_t window_idx = 0;
    int sample_in_window = 0;
    
    for (size_t i = 0; i < nav_samples; ++i) {
        if (sample_in_window >= SAMPLES_PER_WINDOW && window_idx + 1 < eb_per_window.size()) {
            window_idx++;
            sample_in_window = 0;
        }
        
        engine.ins.eb = eb_per_window[window_idx];
        engine.Step_Nav(imu_data[i].wm, imu_data[i].vm);
        sample_in_window++;
        
        double dN = (engine.ins.pos(0) - pos_ref(0)) * glv.Re;
        double dE = (engine.ins.pos(1) - pos_ref(1)) * cos(pos_ref(0)) * glv.Re;
        double drift = sqrt(dN*dN + dE*dE);
        if (drift > max_drift) max_drift = drift;
    }
    
    double nav_hours = nav_samples * ts / 3600.0;
    return max_drift / nav_hours / 1852.0;
}

// =============================================================================
// 主函数
// =============================================================================
int main() {
    cout << "==================================================" << endl;
    cout << "  Strategy Comparison for Backtrack Smoothing" << endl;
    cout << "==================================================" << endl;
    
    // 参数
    string csv_path = "../fog3h.csv";
    double ts = 1.0 / 400.0;
    GLV glv;
    
    Vector3d pos_nj(32.0286 * glv.deg, 118.8533 * glv.deg, 17.0);
    Vector3d att_true = OptimalParams::att_deg * glv.deg;
    Vector3d eb_true = OptimalParams::eb_dph * glv.deg / 3600.0;
    
    Earth eth(glv);
    eth.update(pos_nj, Vector3d::Zero());
    double local_g = eth.gn.norm();
    
    // 加载数据
    cout << "\n[1] Loading data..." << endl;
    auto imu_data = LoadIMUData(csv_path, ts, local_g);
    cout << "    Loaded " << imu_data.size() * ts / 60 << " min" << endl;
    
    // 初始化
    SinsEngine engine(ts);
    engine.Sins_Init(pos_nj, Vector3d::Zero(), Vector3d::Zero(), 
                     Vector3d::Zero(), Vector3d::Zero());
    
    AtomicGyroSimulator atom(att_true, pos_nj, glv);
    
    // 收集窗口数据
    double nav_time = 3600.0;  // 1小时
    double T_WINDOW = 200.0;   // 200秒窗口
    
    cout << "\n[2] Collecting window data (T=" << T_WINDOW << "s)..." << endl;
    auto windows = CollectWindowData(imu_data, atom, ts, nav_time, T_WINDOW, glv);
    cout << "    Windows: " << windows.size() << endl;
    
    // ==========================================================
    // 数据分析
    // ==========================================================
    cout << "\n[3] Data Analysis" << endl;
    cout << "---------------------------------------" << endl;
    cout << "Window | Time(s) |  eb_diff_x  |  eb_diff_y  |  eb_diff_z  | (deg/h)" << endl;
    cout << "---------------------------------------" << endl;
    
    vector<double> t_vec, ex_vec, ey_vec, ez_vec;
    for (size_t w = 0; w < windows.size(); ++w) {
        double t_mid = (windows[w].start_time + windows[w].end_time) / 2;
        t_vec.push_back(t_mid);
        ex_vec.push_back(windows[w].eb_diff(0));
        ey_vec.push_back(windows[w].eb_diff(1));
        ez_vec.push_back(windows[w].eb_diff(2));
        
        Vector3d eb_dph = windows[w].eb_diff / glv.deg * 3600;
        cout << setw(6) << w << " | " 
             << setw(7) << fixed << setprecision(0) << t_mid << " | "
             << setw(11) << setprecision(4) << eb_dph(0) << " | "
             << setw(11) << eb_dph(1) << " | "
             << setw(11) << eb_dph(2) << endl;
    }
    cout << "---------------------------------------" << endl;
    
    // 线性拟合
    double ax, bx, ay, by, az, bz;
    LinearFit(t_vec, ex_vec, ax, bx);
    LinearFit(t_vec, ey_vec, ay, by);
    LinearFit(t_vec, ez_vec, az, bz);
    
    cout << "\nLinear fit: eb(t) = a + b*t" << endl;
    cout << "  X: a=" << ax/glv.deg*3600 << " deg/h, b=" << bx/glv.deg*3600*1000 << "e-3 deg/h/s" << endl;
    cout << "  Y: a=" << ay/glv.deg*3600 << " deg/h, b=" << by/glv.deg*3600*1000 << "e-3 deg/h/s" << endl;
    cout << "  Z: a=" << az/glv.deg*3600 << " deg/h, b=" << bz/glv.deg*3600*1000 << "e-3 deg/h/s" << endl;
    
    // 均值和标准差
    Vector3d mean_eb = Vector3d::Zero();
    for (const auto& w : windows) mean_eb += w.eb_diff;
    mean_eb /= windows.size();
    
    Vector3d std_eb = Vector3d::Zero();
    for (const auto& w : windows) {
        Vector3d d = w.eb_diff - mean_eb;
        std_eb += Vector3d(d(0)*d(0), d(1)*d(1), d(2)*d(2));
    }
    std_eb = (std_eb / windows.size()).cwiseSqrt();
    
    cout << "\nStatistics:" << endl;
    cout << "  Mean eb_diff (deg/h): " << (mean_eb / glv.deg * 3600).transpose() << endl;
    cout << "  Std  eb_diff (deg/h): " << (std_eb / glv.deg * 3600).transpose() << endl;
    cout << "  True eb      (deg/h): " << OptimalParams::eb_dph.transpose() << endl;
    
    // ==========================================================
    // 策略对比
    // ==========================================================
    cout << "\n[4] Strategy Comparison" << endl;
    cout << "---------------------------------------" << endl;
    
    size_t n_windows = windows.size();
    
    // 策略1: 不修正 (eb=0)
    {
        vector<Vector3d> eb_vec(n_windows, Vector3d::Zero());
        double drift = NavigateWithBias(engine, imu_data, att_true, pos_nj, 
                                        eb_vec, ts, T_WINDOW, nav_time, glv);
        cout << "  1. No correction (eb=0):      " << fixed << setprecision(3) << drift << " nm/h" << endl;
    }
    
    // 策略2: 真值零偏
    {
        vector<Vector3d> eb_vec(n_windows, eb_true);
        double drift = NavigateWithBias(engine, imu_data, att_true, pos_nj, 
                                        eb_vec, ts, T_WINDOW, nav_time, glv);
        cout << "  2. True bias (optimal):       " << fixed << setprecision(3) << drift << " nm/h" << endl;
    }
    
    // 策略3: 第一个窗口的 eb_diff
    {
        vector<Vector3d> eb_vec(n_windows, windows[0].eb_diff);
        double drift = NavigateWithBias(engine, imu_data, att_true, pos_nj, 
                                        eb_vec, ts, T_WINDOW, nav_time, glv);
        cout << "  3. First window eb_diff:      " << fixed << setprecision(3) << drift << " nm/h" << endl;
    }
    
    // 策略4: 每个窗口各自的 eb_diff
    {
        vector<Vector3d> eb_vec;
        for (const auto& w : windows) eb_vec.push_back(w.eb_diff);
        double drift = NavigateWithBias(engine, imu_data, att_true, pos_nj, 
                                        eb_vec, ts, T_WINDOW, nav_time, glv);
        cout << "  4. Per-window eb_diff:        " << fixed << setprecision(3) << drift << " nm/h" << endl;
    }
    
    // 策略5: 全局均值
    {
        vector<Vector3d> eb_vec(n_windows, mean_eb);
        double drift = NavigateWithBias(engine, imu_data, att_true, pos_nj, 
                                        eb_vec, ts, T_WINDOW, nav_time, glv);
        cout << "  5. Global mean eb_diff:       " << fixed << setprecision(3) << drift << " nm/h" << endl;
    }
    
    // 策略6: 线性拟合
    {
        vector<Vector3d> eb_vec;
        for (size_t w = 0; w < n_windows; ++w) {
            double t = t_vec[w];
            eb_vec.push_back(Vector3d(ax + bx*t, ay + by*t, az + bz*t));
        }
        double drift = NavigateWithBias(engine, imu_data, att_true, pos_nj, 
                                        eb_vec, ts, T_WINDOW, nav_time, glv);
        cout << "  6. Linear fit eb(t):          " << fixed << setprecision(3) << drift << " nm/h" << endl;
    }
    
    // 策略7: 滑动平均 (前3个窗口)
    {
        vector<Vector3d> eb_vec;
        for (size_t w = 0; w < n_windows; ++w) {
            Vector3d sum = Vector3d::Zero();
            int count = 0;
            for (int k = max(0, (int)w - 2); k <= (int)w; ++k) {
                sum += windows[k].eb_diff;
                count++;
            }
            eb_vec.push_back(sum / count);
        }
        double drift = NavigateWithBias(engine, imu_data, att_true, pos_nj, 
                                        eb_vec, ts, T_WINDOW, nav_time, glv);
        cout << "  7. Sliding avg (N=3):         " << fixed << setprecision(3) << drift << " nm/h" << endl;
    }
    
    // 策略8: 滑动平均 (前5个窗口)
    {
        vector<Vector3d> eb_vec;
        for (size_t w = 0; w < n_windows; ++w) {
            Vector3d sum = Vector3d::Zero();
            int count = 0;
            for (int k = max(0, (int)w - 4); k <= (int)w; ++k) {
                sum += windows[k].eb_diff;
                count++;
            }
            eb_vec.push_back(sum / count);
        }
        double drift = NavigateWithBias(engine, imu_data, att_true, pos_nj, 
                                        eb_vec, ts, T_WINDOW, nav_time, glv);
        cout << "  8. Sliding avg (N=5):         " << fixed << setprecision(3) << drift << " nm/h" << endl;
    }
    
    cout << "---------------------------------------" << endl;
    
    // ==========================================================
    // 保存数据用于绘图
    // ==========================================================
    ofstream fout("eb_drift_analysis.csv");
    fout << "Window,Time_s,eb_x_dph,eb_y_dph,eb_z_dph,fog_z_dph,atom_z_dph" << endl;
    for (size_t w = 0; w < windows.size(); ++w) {
        Vector3d eb_dph = windows[w].eb_diff / glv.deg * 3600;
        fout << w << "," << t_vec[w] << ","
             << eb_dph(0) << "," << eb_dph(1) << "," << eb_dph(2) << ","
             << windows[w].w_fog_mean(2) / glv.deg * 3600 << ","
             << windows[w].w_atom_mean(2) / glv.deg * 3600 << endl;
    }
    fout.close();
    cout << "\nData saved to eb_drift_analysis.csv" << endl;
    
    return 0;
}
