#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include "sins_engine.h"
#include "support.h"

using namespace std;
using namespace Eigen;

int main() {
    string csv_path = "../fog3h.csv"; 
    double ts = 1.0 / 400.0; 
    GLV glv;

    // ---------------------------------------------------------
    // [CRITICAL FIX] 计算初始位置的"当地重力"，用于归一化
    // ---------------------------------------------------------
    // 初始位置：南京
    Vector3d pos_nj(32.0286 * glv.deg, 118.8533 * glv.deg, 17.0); 
    
    // 借用 Earth 类算一下当地重力
    Earth eth_temp(glv);
    eth_temp.update(pos_nj, Vector3d::Zero());
    double local_g = eth_temp.gn.norm(); // 理论值约为 9.7932
    
    cout << "[Config] Standard Gravity (g0): " << glv.g0 << endl;
    cout << "[Config] Local Gravity (gn):    " << local_g << " (Difference: " << (local_g - glv.g0) << ")" << endl;
    
    // ---------------------------------------------------------
    // 配置与加载
    // ---------------------------------------------------------
    // 先跑 600s 快速验证。确认漂移率下降后，可改为 10800.0 跑全量
    double t_coarse = 60.0;
    double t_fine = 600.0;  
    
    // [关键] 传入 local_g 进行归一化
    auto sim = LoadAndSplitData(csv_path, ts, local_g, t_coarse, t_fine);
    
    if (sim.data_fine.empty()) { cerr << "Data load failed or too short." << endl; return -1; }
    
    SinsEngine engine(ts);
    // 初始化 INS 状态
    engine.Sins_Init(pos_nj, Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero());

    // ---------------------------------------------------------
    // 1. 粗对准
    // ---------------------------------------------------------
    cout << ">>> [1/3] Running Coarse Alignment (" << sim.time_coarse << "s)..." << endl;
    //engine.Run_Coarse_Phase(sim.data_coarse);
    Vector3d att_truth(0.03404, 0.31606, 0.61788); // 你的真值
    engine.res_coarse.valid = true;
    engine.res_coarse.att = att_truth * glv.deg; // 注意单位
    engine.res_coarse.vel.setZero();
    engine.res_coarse.pos = pos_nj;
    engine.res_coarse.eb.setZero();
    engine.res_coarse.db.setZero();
    if (engine.res_coarse.valid) {
        Vector3d att = engine.res_coarse.att * glv.rad;
        cout << "    Result (deg): P=" << att(0) << ", R=" << att(1) << ", Y=" << att(2) << endl;
    } else {
        cerr << "    Coarse Alignment Failed!" << endl; return -1;
    }

    // ---------------------------------------------------------
    // 2. 精对准 (Fine Alignment) - 开环估计
    // ---------------------------------------------------------
    cout << ">>> [2/3] Running Fine Alignment (" << sim.time_fine << "s)..." << endl;
    
    ofstream f_log("align_log.txt");
    f_log << "Time,Att0,Att1,Att2,Vel0,Vel1,Vel2,Eb0,Eb1,Eb2,Db0,Db1,Db2" << endl;

    {
        engine.res_fine.valid = false;
        AlignResult* start = &engine.res_coarse;
        
        // 重置 INS 状态
        engine.eth.update(start->pos, start->vel);
        engine.ins = INSState(start->att, start->vel, start->pos, engine.ins.ts, engine.eth);
        engine.ins.set_bias(start->eb, start->db);
        
        // 初始化 KF
        engine.kf.Init(engine.ins.ts, glv, 
            engine.kf_cfg.phi_init_err(0), engine.kf_cfg.web_psd, engine.kf_cfg.wdb_psd, 
            engine.kf_cfg.eb_sigma, engine.kf_cfg.db_sigma, engine.kf_cfg.wvn_err(0));
            
        int progress_step = sim.data_fine.size() / 10;
        for (size_t i = 0; i < sim.data_fine.size(); ++i) {
            engine.Step_Fine(sim.data_fine[i].wm, sim.data_fine[i].vm);
            
            // 降频记录日志
            if (i % 400 == 0) { 
                Vector3d att = engine.GetAttDeg();
                Vector3d vel = engine.GetVel();
                Vector3d eb = engine.GetBiasGyro(); 
                Vector3d db = engine.GetBiasAcc(); 
                f_log << sim.data_fine[i].t << "," 
                      << att(0) << "," << att(1) << "," << att(2) << ","
                      << vel(0) << "," << vel(1) << "," << vel(2) << ","
                      << eb(0) << "," << eb(1) << "," << eb(2) << ","
                      << db(0) << "," << db(1) << "," << db(2) << endl;
            }
            if (i > 0 && i % progress_step == 0) cout << "    Progress: " << (i*100/sim.data_fine.size()) << "%" << endl;
        }
        engine.Finish_Fine(); 
    }
    f_log.close();

    // ---------------------------------------------------------
    // 3. 导航验证 (使用剩余数据)
    // ---------------------------------------------------------
    size_t start_idx = static_cast<size_t>(t_fine / ts); 
    
    if (start_idx < sim.data_nav.size()) {
        double nav_duration = 3600.0; // 跑1小时导航 (数据不够则自动截止)
        cout << ">>> [3/3] Running Nav Verification (" << nav_duration << "s limit)..." << endl;
        
        ofstream f_nav("nav_log.txt");
        f_nav << "Time,Lat,Lon,H,Ve,Vn,Vu" << endl;

        // 初始化导航状态
        engine.Run_Nav_Phase(std::vector<IMUData>()); 
        
        Vector3d p0 = engine.GetPosDeg();
        
        size_t end_idx = start_idx + static_cast<size_t>(nav_duration / ts);
        if (end_idx > sim.data_nav.size()) end_idx = sim.data_nav.size();

        int log_step = 0;
        for (size_t i = start_idx; i < end_idx; ++i) {
            const auto& d = sim.data_nav[i];
            engine.Step_Nav(d.wm, d.vm); 
            
            if (++log_step >= 400) { 
                log_step = 0;
                Vector3d p = engine.GetPosDeg();
                Vector3d v = engine.GetVel();
                f_nav << d.t << "," << p(0) << "," << p(1) << "," << p(2) << "," 
                      << v(0) << "," << v(1) << "," << v(2) << endl;
            }
        }
        f_nav.close();

        Vector3d p1 = engine.GetPosDeg();
        double drift = (p1.head<2>() - p0.head<2>()).norm() * 111320.0; 
        
        cout << "------------------------------------------------" << endl;
        cout << "Init Pos: " << p0.transpose() << endl;
        cout << "End Pos:  " << p1.transpose() << endl;
        cout << "Total Horizontal Drift: " << drift << " m" << endl;
        
        double actual_time_h = (end_idx - start_idx) * ts / 3600.0;
        if(actual_time_h > 0)
             cout << "Drift Rate: " << (drift / 1852.0) / actual_time_h << " nm/h" << endl;
    }

    return 0;
}