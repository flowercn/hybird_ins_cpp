#ifndef SUPPORT_H
#define SUPPORT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <Eigen/Dense>
#include "../psins/sins_engine.h"  

using namespace std;
using namespace Eigen;

// --- [原有功能] 数据加载 ---
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

// --- [原有功能] 坐标转换 ---
Vector3d LocalToGeo(const Vector3d& xyz_enu, const Vector3d& pos_ref_rad, const GLV& glv) {
    double lat = pos_ref_rad(0); double h = pos_ref_rad(2);
    double sl = sin(lat); double cl = cos(lat); 
    double e2 = glv.e2; double Re = glv.Re;
    double sq = sqrt(1.0 - e2 * sl * sl);
    double RN = Re / sq; double RM = RN * (1.0 - e2) / (1.0 - e2 * sl * sl); 
    return Vector3d(lat + xyz_enu(1)/(RM+h), pos_ref_rad(1) + xyz_enu(0)/((RN+h)*cl), h + xyz_enu(2));
}

// --- [新增功能] 惯导引擎一键初始化 ---
// 将 main 函数中繁琐的配置代码封装在此
SinsEngine CreateConfiguredSins(double ts, const GLV& glv, double t_coarse, double t_fine) {
    AlignConfig align_cfg;
    
    // 姿态参考 (Attitude Reference)
    align_cfg.att_ref << -0.03258, 0.20927, 0.62977; 
    align_cfg.att_ref *= glv.deg;  
    
    // 位置参考 (Position Reference)
    align_cfg.pos_ref << 32.0286 * glv.deg, 118.8533 * glv.deg, 17.0; 
    
    // 初始误差设置 (Initial Errors)
    align_cfg.phi_init_err << 0.01, 0.01, 0.1;
    align_cfg.phi_init_err *= glv.deg; 
    align_cfg.wvn_err << 0.001, 0.001, 0.001;
    
    // IMU 噪声参数 (Sensor Noise)
    align_cfg.eb_sigma = 0.2 * glv.dph;          
    align_cfg.db_sigma = 100 * glv.ug;           
    align_cfg.web_psd  = 0.021 * glv.dpsh;       
    align_cfg.wdb_psd  = 10 * glv.ugpsHz;

    // 实例化并初始化引擎
    SinsEngine sins(ts);
    sins.Init(align_cfg, t_coarse, t_fine);
    
    return sins; // 返回初始化好的对象
}

#endif // SUPPORT_H