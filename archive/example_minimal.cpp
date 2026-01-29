/**
 * @file example_minimal.cpp
 * @brief 最简化的纯惯导示例 - 不到 50 行代码
 * 
 * 展示如何用最少的代码完成：
 *   1. 加载数据
 *   2. 混合对准 (KF姿态 + 几何均值零偏)
 *   3. 纯惯导导航
 */

#include <iostream>
#include <iomanip>
#include "sins_engine.h"
#include "support.h"

using namespace std;
using namespace Eigen;

int main() {
    // ========== 参数 ==========
    double ts = 1.0 / 400.0;
    GLV glv;
    Vector3d pos(32.0286 * glv.deg, 118.8533 * glv.deg, 17.0);  // 南京
    
    // 当地重力
    Earth eth(glv);
    eth.update(pos, Vector3d::Zero());
    
    // ========== 加载数据 ==========
    auto imu_data = LoadIMUData("../fog3h.csv", ts, eth.gn.norm());
    cout << "Loaded " << imu_data.size() * ts / 60 << " min data" << endl;
    
    // ========== 初始化 + 对准 ==========
    SinsEngine ins(ts);
    ins.Sins_Init(pos, Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero());
    
    // 一键对准：60s粗对准 + 1h精对准
    HybridAlignConfig cfg;
    cfg.t_coarse = 60;
    cfg.t_fine = 3600;
    cfg.verbose = false;  // 静默模式
    
    auto result = ins.Run_HybridAlign(imu_data, cfg);
    if (!result.valid) { cerr << "Align failed!" << endl; return -1; }
    
    cout << "Alignment done. Attitude (deg): " << (result.att / glv.deg).transpose() << endl;
    cout << "Bias (deg/h): " << (result.eb / glv.deg * 3600).transpose() << endl;
    
    // ========== 导航 ==========
    // 重置到对准后状态
    ins.eth.update(pos, Vector3d::Zero());
    ins.ins = INSState(result.att, Vector3d::Zero(), pos, ts, ins.eth);
    ins.ins.set_bias(result.eb, result.db);
    
    // 从对准结束处开始导航
    size_t nav_start = static_cast<size_t>(cfg.t_fine / ts);
    double max_drift = 0;
    
    for (size_t i = nav_start; i < imu_data.size(); ++i) {
        ins.Step_Nav(imu_data[i].wm, imu_data[i].vm);
        
        double dN = (ins.ins.pos(0) - pos(0)) * glv.Re;
        double dE = (ins.ins.pos(1) - pos(1)) * cos(pos(0)) * glv.Re;
        double drift = sqrt(dN*dN + dE*dE);
        if (drift > max_drift) max_drift = drift;
    }
    
    double nav_time = (imu_data.size() - nav_start) * ts;
    cout << "\nNav " << nav_time/60 << " min, max drift: " << max_drift << " m ("
         << fixed << setprecision(3) << max_drift / (nav_time/3600) / 1852 << " nm/h)" << endl;
    
    return 0;
}
