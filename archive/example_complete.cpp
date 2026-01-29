/**
 * @file example_complete.cpp
 * @brief 完整的纯惯导示例 - 展示所有接口用法
 * 
 * 功能演示：
 *   1. 数据加载与配置
 *   2. 混合对准 (详细输出)
 *   3. 纯惯导导航 (实时输出)
 *   4. 结果分析与对比
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include "sins_engine.h"
#include "support.h"

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[]) {
    // ========== 1. 参数配置 ==========
    string csv_path = "../fog3h.csv";
    double ts = 1.0 / 400.0;
    
    // 命令行参数
    double t_coarse = 60.0;       // 粗对准时间 (s)
    double t_fine = 3600.0;       // 精对准时间 (s)
    double t_nav = 3600.0;        // 导航时间 (s)
    double eb_sigma = 0.003;      // Allan 零偏稳定性 (deg/h)
    
    if (argc > 1) t_fine = atof(argv[1]);  // 可指定精对准时间
    if (argc > 2) t_nav = atof(argv[2]);   // 可指定导航时间
    
    GLV glv;
    Vector3d pos_init(32.0286 * glv.deg, 118.8533 * glv.deg, 17.0);
    
    cout << "================================================" << endl;
    cout << "    Pure INS Navigation Demo" << endl;
    cout << "================================================" << endl;
    cout << "Position: Nanjing (32.03°N, 118.85°E)" << endl;
    cout << "Sample rate: " << 1.0/ts << " Hz" << endl;
    cout << "Coarse align: " << t_coarse << " s" << endl;
    cout << "Fine align:   " << t_fine/60 << " min" << endl;
    cout << "Navigation:   " << t_nav/60 << " min" << endl;
    cout << "================================================\n" << endl;
    
    // ========== 2. 加载数据 ==========
    cout << "[1] Loading IMU data..." << endl;
    
    Earth eth_temp(glv);
    eth_temp.update(pos_init, Vector3d::Zero());
    double local_g = eth_temp.gn.norm();
    
    auto imu_data = LoadIMUData(csv_path, ts, local_g);
    if (imu_data.empty()) {
        cerr << "Failed to load data!" << endl;
        return -1;
    }
    
    double total_time = imu_data.size() * ts;
    cout << "    Loaded " << imu_data.size() << " samples (" << total_time/60 << " min)" << endl;
    cout << "    Local gravity: " << local_g << " m/s^2" << endl;
    
    // ========== 3. 初始化惯导 ==========
    cout << "\n[2] Initializing INS engine..." << endl;
    
    SinsEngine engine(ts);
    engine.Sins_Init(pos_init, Vector3d::Zero(), Vector3d::Zero(), 
                     Vector3d::Zero(), Vector3d::Zero());
    
    // ========== 4. 混合对准 ==========
    cout << "\n[3] Running Hybrid Alignment..." << endl;
    
    HybridAlignConfig align_cfg;
    align_cfg.t_coarse = t_coarse;
    align_cfg.t_fine = t_fine;
    align_cfg.eb_sigma_allan = eb_sigma;
    align_cfg.verbose = true;
    
    auto align = engine.Run_HybridAlign(imu_data, align_cfg);
    
    if (!align.valid) {
        cerr << "Alignment failed!" << endl;
        return -1;
    }
    
    // ========== 5. 纯惯导导航 ==========
    cout << "\n[4] Starting Pure INS Navigation..." << endl;
    
    // 重置 INS 状态
    engine.eth.update(pos_init, Vector3d::Zero());
    engine.ins = INSState(align.att, Vector3d::Zero(), pos_init, ts, engine.eth);
    engine.ins.set_bias(align.eb, align.db);
    
    // 导航数据起点
    size_t nav_start = static_cast<size_t>(t_fine / ts);
    size_t nav_end = min(nav_start + static_cast<size_t>(t_nav / ts), imu_data.size());
    size_t nav_count = nav_end - nav_start;
    
    // 轨迹记录
    vector<double> time_log, drift_log;
    double max_drift = 0;
    
    cout << "\n    Time(min)  Drift(m)  Pos(Lat,Lon)" << endl;
    cout << "    -----------------------------------------" << endl;
    
    int print_interval = nav_count / 10;
    
    for (size_t i = nav_start; i < nav_end; ++i) {
        engine.Step_Nav(imu_data[i].wm, imu_data[i].vm);
        
        // 计算漂移
        double dN = (engine.ins.pos(0) - pos_init(0)) * glv.Re;
        double dE = (engine.ins.pos(1) - pos_init(1)) * cos(pos_init(0)) * glv.Re;
        double drift = sqrt(dN*dN + dE*dE);
        
        if (drift > max_drift) max_drift = drift;
        
        // 记录
        size_t idx = i - nav_start;
        if (idx % 1000 == 0) {
            time_log.push_back((idx * ts) / 60.0);
            drift_log.push_back(drift);
        }
        
        // 输出
        if (print_interval > 0 && idx > 0 && idx % print_interval == 0) {
            double t_min = (idx * ts) / 60.0;
            Vector3d pos_deg = engine.GetPosDeg();
            cout << "    " << fixed << setw(8) << setprecision(1) << t_min
                 << setw(10) << setprecision(1) << drift
                 << "    (" << setprecision(5) << pos_deg(0) << ", " << pos_deg(1) << ")" << endl;
        }
    }
    
    // ========== 6. 结果汇总 ==========
    double nav_time_actual = nav_count * ts;
    double drift_rate = max_drift / (nav_time_actual / 3600.0) / 1852.0;
    
    cout << "\n================================================" << endl;
    cout << "    Navigation Summary" << endl;
    cout << "================================================" << endl;
    cout << "Navigation time:  " << nav_time_actual/60 << " min" << endl;
    cout << "Max drift:        " << fixed << setprecision(1) << max_drift << " m" << endl;
    cout << "Drift rate:       " << setprecision(3) << drift_rate << " nm/h" << endl;
    cout << endl;
    cout << "Final attitude (deg):" << endl;
    cout << "  Pitch: " << engine.GetAttDeg()(0) << endl;
    cout << "  Roll:  " << engine.GetAttDeg()(1) << endl;
    cout << "  Yaw:   " << engine.GetAttDeg()(2) << endl;
    cout << endl;
    cout << "Final position:" << endl;
    Vector3d pos_final = engine.GetPosDeg();
    cout << "  Lat: " << pos_final(0) << "°" << endl;
    cout << "  Lon: " << pos_final(1) << "°" << endl;
    cout << "================================================" << endl;
    
    // ========== 7. 保存轨迹 (可选) ==========
    ofstream fout("nav_trajectory.csv");
    if (fout.is_open()) {
        fout << "time_min,drift_m" << endl;
        for (size_t i = 0; i < time_log.size(); ++i) {
            fout << time_log[i] << "," << drift_log[i] << endl;
        }
        fout.close();
        cout << "\nTrajectory saved to nav_trajectory.csv" << endl;
    }
    
    return 0;
}
