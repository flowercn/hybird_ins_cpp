#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <memory>
#include <algorithm> 
#include <Eigen/Dense>

#include "../psins/sins_engine.h"  
#include "cai_sim.h"
#include "fgo_core.h"

using namespace std;
using namespace Eigen;

std::vector<IMUData> load_and_process_data(const std::string& filename, double ts, double target_g);
Vector3d LocalToGeo(const Vector3d& xyz_enu, const Vector3d& pos_ref_rad, const GLV& glv);

int main() {
    GLV glv; 
    string fog_path = "../fog.csv"; 
    double ts = 1.0 / 400.0;

    // [Step 0] 严格物理缩放：确保数据重力 = 算法重力
    cout << "[DEBUG] Loading and scaling data..." << endl;
    auto all_data = load_and_process_data(fog_path, ts, glv.g0);
    if (all_data.empty()) return -1;

    // SINS 粗对准参数
    AlignConfig align_cfg;
    align_cfg.att_ref << -0.03258, 0.20927, 0.62977; align_cfg.att_ref *= glv.deg;  
    align_cfg.pos_ref << 32.0286 * glv.deg, 118.8533 * glv.deg, 17.0; 
    align_cfg.phi_init_err << 0.01, 0.01, 0.1; align_cfg.phi_init_err *= glv.deg; 
    align_cfg.wvn_err << 0.001, 0.001, 0.001;
    align_cfg.eb_sigma = 0.2 * glv.dph; align_cfg.db_sigma = 100 * glv.ug;          
    align_cfg.web_psd  = 0.021 * glv.dpsh; align_cfg.wdb_psd  = 10 * glv.ugpsHz;

    SinsEngine sins(ts);
    double align_duration = 300.0; // 对准 300秒
    sins.Init(align_cfg, align_duration);

    CAIGSimulator cai;
    cai.Init(align_duration);

    // FGO 配置
    FGOSettings fgo_opts;
    fgo_opts.fog_dt = ts;
    fgo_opts.node_interval = 1.0; 
    fgo_opts.gravity = glv.g0; // 物理模型一致
    
    // 设置噪声参数
    fgo_opts.gyro_noise = 0.1 * glv.deg;      
    fgo_opts.gyro_bias_noise = 1.0e-5;         
    fgo_opts.accel_bias_noise = 1.0e-8; // 锁死加计Bias
    fgo_opts.cai_rot_noise = 2e-4 * glv.deg * sqrt(2.0); 
    
    FGOEngine fgo(fgo_opts);
    bool fgo_started = false;

    // 日志
    std::ofstream fgo_log("res_fgo_final.csv");
    fgo_log << "time,lat,lon,h,vn,ve,vd,roll,pitch,yaw,bg_x,bg_y,bg_z,ba_x,ba_y,ba_z" << endl;

    cout << "🚀 Starting Relay Simulation (SINS Align -> FGO Nav)..." << endl;
    
    for (const auto& d : all_data) {
        // SINS 始终在运行，用于提供对准结果和作为 CAIG 的“真值驱动”
        sins.Step(d.wm, d.vm, d.t);

        if (d.t >= align_duration) {
            
            if (!fgo_started) {
                // [Step 1: 完美继承]
                // 1. 姿态：用 SINS 算好的精对准姿态
                Vector3d att_init = sins.GetAttDeg() * glv.deg; 
                // 2. 速度/位置：按你要求置零 (本地导航系原点)
                Vector3d vel_init = Vector3d::Zero(); 
                Vector3d pos_init = Vector3d::Zero(); 
                // 3. 零偏：用 SINS 估计出的零偏 (eb, db)
                Vector3d bg_init = sins.GetBiasGyro();
                Vector3d ba_init = sins.GetBiasAcc();

                // 初始化 FGO
                fgo.initialize(d.t, att_init, vel_init, pos_init, bg_init, ba_init);
                fgo_started = true;
            }

            // 原子陀螺仿真
            Vector3d cai_omega;
            bool cai_has_data = cai.Update(d.t, sins.GetQnb(), cai_omega);

            Eigen::Vector3d gyro_input = Eigen::Vector3d::Zero(); // 平台稳定
            Eigen::Vector3d acc_input = d.vm / ts; 

            // [Step 2: 严禁 ZUPT]
            // force_align = false -> 纯积分 + 原子观测
            bool new_node = fgo.process_imu(gyro_input, acc_input, ts, false, false);
            
            if (cai_has_data) {
                fgo.add_cai_measurement(d.t, cai_omega);
                // 进度打印
                static double last_print = 0;
                if (d.t - last_print > 10.0) {
                    cout << "\r[Nav] t=" << fixed << setprecision(1) << d.t << "s" << flush;
                    last_print = d.t;
                }
            }

            if (new_node) {
                NavStateResult res = fgo.get_result();
                Vector3d ypr = res.pose.rotation().ypr(); 
                Vector3d att_deg(ypr(1), ypr(2), ypr(0)); att_deg *= glv.rad;
                
                // 将 FGO 的局部位置加到参考点上
                Vector3d local_xyz(res.pose.x(), res.pose.y(), res.pose.z());
                Vector3d geo_pos = LocalToGeo(local_xyz, align_cfg.pos_ref, glv);
                
                Vector3d bg = res.bias.gyroscope() * glv.rad * 3600.0; 
                Vector3d ba = res.bias.accelerometer() / glv.ug;       

                fgo_log << fixed << setprecision(9) 
                        << res.time << ","
                        << geo_pos(0)*glv.rad << "," << geo_pos(1)*glv.rad << "," << geo_pos(2) << ","
                        << res.vel(0) << "," << res.vel(1) << "," << res.vel(2) << ","
                        << att_deg(0) << "," << att_deg(1) << "," << att_deg(2) << ","
                        << bg(0) << "," << bg(1) << "," << bg(2) << ","
                        << ba(0) << "," << ba(1) << "," << ba(2) << endl;
            }
        }
    }

    cout << "\n✅ Simulation Finished." << endl;
    return 0;
}

// -----------------------------------------------------------------------
// 保持原样的数据加载函数
// -----------------------------------------------------------------------
std::vector<IMUData> load_and_process_data(const std::string& filename, double ts, double target_g) {
    std::vector<Vector3d> raw_gyro, raw_acc;
    std::ifstream file(filename);
    std::string line;
    
    if (!file.is_open()) { cerr << "Error: Cannot open " << filename << endl; return {}; }
    
    while (std::getline(file, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string cell;
        std::vector<double> vals;
        try {
            while (std::getline(ss, cell, ',')) {
                cell.erase(std::remove(cell.begin(), cell.end(), '\r'), cell.end()); 
                if (!cell.empty()) vals.push_back(std::stod(cell));
            }
        } catch (...) { continue; } 
        
        if (vals.size() >= 6) {
            raw_gyro.push_back(Vector3d(vals[0], vals[1], vals[2]));
            raw_acc.push_back(Vector3d(vals[3], vals[4], vals[5]));
        }
    }

    if (raw_gyro.empty()) return {};

    Vector3d acc_sum = Vector3d::Zero();
    for (const auto& a : raw_acc) acc_sum += a;
    Vector3d acc_mean_raw = acc_sum / raw_acc.size();
    
    double norm_meas = acc_mean_raw.norm();
    double scale_ratio = (norm_meas > 1e-5) ? target_g / norm_meas : 1.0; 
    
    cout << "  Scale Ratio: " << scale_ratio << endl;

    std::vector<IMUData> data;
    double t_curr = 0.0;
    for (size_t i = 0; i < raw_gyro.size(); ++i) {
        IMUData d;
        d.wm = raw_gyro[i] * ts;          
        d.vm = (raw_acc[i] * scale_ratio) * ts; 
        d.t = t_curr;
        data.push_back(d);
        t_curr += ts;
    }
    return data;
}

Vector3d LocalToGeo(const Vector3d& xyz_enu, const Vector3d& pos_ref_rad, const GLV& glv) {
    double lat = pos_ref_rad(0); double h = pos_ref_rad(2);
    double sl = sin(lat); double cl = cos(lat); 
    double e2 = glv.e2; double Re = glv.Re;
    double sq = sqrt(1.0 - e2 * sl * sl);
    double RN = Re / sq; double RM = RN * (1.0 - e2) / (1.0 - e2 * sl * sl); 
    return Vector3d(lat + xyz_enu(1)/(RM+h), pos_ref_rad(1) + xyz_enu(0)/((RN+h)*cl), h + xyz_enu(2));
}