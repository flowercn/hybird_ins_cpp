/**
 * @file theory_test.cpp
 * @brief 简单理论验证 - 直接用误差公式计算
 * 
 * 验证思路：
 *   1. 从真实数据提取噪声特性
 *   2. 用公式计算理论极限
 *   3. 对比实际导航结果
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <random>
#include <numeric>

#include "sins_engine.h"
#include "support.h"

using namespace std;
using namespace Eigen;

int main() {
    cout << "============================================================\n";
    cout << "  惯导理论极限验证\n";
    cout << "============================================================\n";
    
    GLV glv;
    double ts = 1.0 / 400.0;
    double T = 3600.0;  // 1小时
    double g = 9.8;
    
    // ========== 1. 加载真实数据 ==========
    cout << "\n[1] 加载真实 FOG 数据...\n";
    
    Vector3d pos0(32.0286 * glv.deg, 118.8533 * glv.deg, 17.0);
    Earth eth(glv);
    eth.update(pos0, Vector3d::Zero());
    double local_g = eth.gn.norm();
    
    auto imu_data = LoadIMUData("../fog3h.csv", ts, local_g);
    cout << "    样本数: " << imu_data.size() << "\n";
    
    // 刻度因数修正
    double K = 1.000372;
    for (auto& imu : imu_data) imu.wm *= K;
    
    // ========== 2. 提取陀螺噪声特性 ==========
    cout << "\n[2] 分析陀螺噪声特性...\n";
    
    // 计算陀螺数据的标准差（作为噪声估计）
    size_t n_align = 2 * 3600 * 400;  // 2小时对准数据
    
    // 先计算均值
    Vector3d w_mean = Vector3d::Zero();
    for (size_t i = 0; i < n_align; ++i) {
        w_mean += imu_data[i].wm / ts;  // 转为角速度
    }
    w_mean /= n_align;
    
    // 计算方差
    Vector3d w_var = Vector3d::Zero();
    for (size_t i = 0; i < n_align; ++i) {
        Vector3d w = imu_data[i].wm / ts;
        Vector3d diff = w - w_mean;
        w_var += diff.cwiseProduct(diff);
    }
    w_var /= n_align;
    Vector3d w_std = w_var.cwiseSqrt();
    
    // 转换为 ARW (deg/sqrt(h))
    // σ_w (rad/s) -> ARW = σ_w * sqrt(fs) * 180/π * sqrt(3600)
    double fs = 1.0 / ts;
    Vector3d arw = w_std * sqrt(fs) * 180.0 / M_PI * sqrt(3600.0);
    
    cout << "    陀螺噪声 σ_w (deg/h): " 
         << w_std(0) * 180/M_PI * 3600 << ", "
         << w_std(1) * 180/M_PI * 3600 << ", "
         << w_std(2) * 180/M_PI * 3600 << "\n";
    cout << "    等效 ARW (deg/√h): "
         << arw(0) << ", " << arw(1) << ", " << arw(2) << "\n";
    
    // ========== 3. 理论公式计算 ==========
    cout << "\n[3] 理论公式计算...\n";
    
    // 使用平均 ARW
    double arw_avg = (arw(0) + arw(1) + arw(2)) / 3.0;
    double arw_rad = arw_avg * M_PI / 180.0;  // deg/√h -> rad/√h
    
    // 姿态误差 σ_θ(T) = ARW * √T (T in hours)
    double sigma_theta = arw_rad * sqrt(T / 3600.0);  // rad
    
    // 位置误差（简化模型，系数 1/6 来自随机积分）
    double drift_theory = (1.0/6.0) * g * sigma_theta * T * T;  // m
    
    cout << "    平均 ARW: " << arw_avg << " deg/√h\n";
    cout << "    1h 姿态误差 σ_θ: " << sigma_theta * 180/M_PI << " deg\n";
    cout << "    理论位置漂移: " << drift_theory << " m (" 
         << drift_theory / 1852.0 << " nm)\n";
    
    // ========== 4. 传统 INS 实际结果 ==========
    cout << "\n[4] 传统 INS 导航...\n";
    
    // 使用之前搜索的最优参数
    Vector3d att_deg(0.0344, 0.330, 0.590);
    Vector3d att0 = att_deg * glv.deg;
    Vector3d eb_dph(0, 0, 0.019);
    Vector3d eb = eb_dph * glv.dph;
    
    SinsEngine ins(ts);
    ins.eth.update(pos0, Vector3d::Zero());
    ins.ins = INSState(att0, Vector3d::Zero(), pos0, ts, ins.eth);
    ins.ins.set_bias(eb, Vector3d::Zero());
    
    // 先处理对准时段的数据（让状态演化）
    cout << "    处理对准时段 (" << n_align << " samples)...\n";
    for (size_t i = 0; i < n_align; ++i) {
        ins.Step_Nav(imu_data[i].wm, imu_data[i].vm);
    }
    
    // 记录对准结束时的位置作为基准
    Vector3d pos_nav_start = ins.ins.pos;
    cout << "    对准后位置: " << pos_nav_start(0)/glv.deg << ", " 
         << pos_nav_start(1)/glv.deg << ", " << pos_nav_start(2) << "\n";
    
    // 导航
    cout << "    开始导航...\n";
    size_t nav_start = n_align;
    double max_drift = 0;
    
    for (size_t i = nav_start; i < imu_data.size(); ++i) {
        ins.Step_Nav(imu_data[i].wm, imu_data[i].vm);
        
        double dN = (ins.ins.pos(0) - pos_nav_start(0)) * glv.Re;
        double dE = (ins.ins.pos(1) - pos_nav_start(1)) * cos(pos_nav_start(0)) * glv.Re;
        double drift = sqrt(dN*dN + dE*dE);
        if (drift > max_drift) max_drift = drift;
        
        // Debug: 打印几个时刻的漂移
        if (i == nav_start + 60*400 || i == nav_start + 30*60*400 || i == nav_start + 59*60*400) {
            cout << "    t=" << (i-nav_start)*ts/60 << " min, drift=" << drift << " m\n";
        }
    }
    
    double nav_time = (imu_data.size() - nav_start) * ts;
    double drift_rate = max_drift / (nav_time / 3600.0) / 1852.0;
    
    cout << "    导航时间: " << nav_time / 60 << " min\n";
    cout << "    最大漂移: " << max_drift << " m\n";
    cout << "    漂移率: " << drift_rate << " nm/h\n";
    
    // ========== 5. 对比分析 ==========
    cout << "\n[5] 对比分析\n";
    cout << "    ============================================\n";
    cout << "    理论 (纯ARW):        " << fixed << setprecision(4) 
         << drift_theory / 1852.0 << " nm/h\n";
    cout << "    实际 (传统INS):      " << drift_rate << " nm/h\n";
    cout << "    实际/理论 = " << setprecision(2) 
         << drift_rate / (drift_theory / 1852.0) << "x\n";
    cout << "    ============================================\n";
    
    // ========== 6. 反推误差源 ==========
    cout << "\n[6] 反推主导误差源\n";
    
    // 如果主要是 ARW，漂移 ∝ T^2.5
    // 如果主要是常值姿态误差，漂移 ∝ T^2
    // 如果主要是零偏，漂移 ∝ T^3
    
    // 从实际漂移反推等效常值姿态误差
    double theta_equiv = max_drift / (0.5 * g * T * T);
    cout << "    等效常值姿态误差: " << theta_equiv * 180/M_PI * 3600 << " arcsec\n";
    cout << "                      = " << theta_equiv * 180/M_PI << " deg\n";
    
    // 从实际漂移反推等效零偏
    double eb_equiv = theta_equiv / T;  // rad/s
    cout << "    等效零偏误差: " << eb_equiv * 180/M_PI * 3600 << " deg/h\n";
    
    // ========== 7. 结论 ==========
    cout << "\n[7] 结论\n";
    cout << "    --------------------------------------------\n";
    
    if (drift_rate > drift_theory / 1852.0 * 1.5) {
        cout << "    主导误差源: 常值误差 (姿态/零偏估计)\n";
        cout << "    改进方向: 更精确的对准和零偏估计\n";
    } else {
        cout << "    主导误差源: 随机噪声 (ARW)\n";
        cout << "    已接近理论极限！\n";
    }
    cout << "    --------------------------------------------\n";
    
    // ========== 8. 93米溯源 ==========
    cout << "\n[8] 关于 93 米的说明\n";
    
    // 93米对应什么 ARW?
    double drift_93 = 93;
    // drift = (1/6) * g * ARW * T^2.5 / sqrt(3600)
    // ARW = drift * 6 * sqrt(3600) / (g * T^2.5)
    double T_h = T / 3600.0;
    double arw_for_93 = drift_93 * 6.0 / (g * pow(T, 2) * sqrt(T_h)) * 180/M_PI;
    
    cout << "    93m 对应 ARW ≈ " << arw_for_93 << " deg/√h\n";
    cout << "    实际 ARW ≈ " << arw_avg << " deg/√h\n";
    cout << "    差距约 " << arw_avg / arw_for_93 << " 倍\n";
    
    cout << "\n============================================================\n";
    
    return 0;
}
