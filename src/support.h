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

// ========== 简化加载接口 ==========
// 只加载原始 IMU 数据，不做切片（用于 HybridAlign 自动切片）
inline std::vector<IMUData> LoadIMUData(const std::string& filename, double ts, double target_g) {
    std::vector<IMUData> data;
    std::vector<Vector3d> raw_gyro, raw_acc;
    
    std::ifstream file(filename);
    if (!file.is_open()) { 
        cerr << "Error: Cannot open " << filename << endl; 
        return data; 
    }
    
    std::string line;
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
    
    if (raw_gyro.empty()) return data;
    
    // 加速度计归一化
    Vector3d acc_sum = Vector3d::Zero();
    for (const auto& a : raw_acc) acc_sum += a;
    double norm_meas = (acc_sum / raw_acc.size()).norm();
    double scale = (norm_meas > 1e-5) ? target_g / norm_meas : 1.0;
    
    // 构建数据
    double t = 0;
    for (size_t i = 0; i < raw_gyro.size(); ++i) {
        IMUData d;
        d.wm = raw_gyro[i] * ts;
        d.vm = raw_acc[i] * scale * ts;
        d.t = t;
        data.push_back(d);
        t += ts;
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

#endif // SUPPORT_H


// struct SimulationData {
//     // 粗对准数据包
//     std::vector<IMUData> data_coarse;
//     double time_coarse; // 秒

//     // 精对准数据包
//     std::vector<IMUData> data_fine;
//     double time_fine;   // 秒

//     // 导航数据包 (通常是全量数据)
//     std::vector<IMUData> data_nav;
//     double time_nav;    // 秒
// };

// SimulationData LoadAndSplitData(const std::string& filename, double ts, double target_g, 
//                                 double req_t_coarse, double req_t_fine) 
// {
//     SimulationData sim_data;
//     std::vector<Vector3d> raw_gyro, raw_acc;
    
//     // 1. 读取 CSV 文件 复制前六列到raw_gyro和raw_acc
//     std::ifstream file(filename);
//     std::string line;
//     if (!file.is_open()) { cerr << "Error: Cannot open " << filename << endl; return sim_data; }
    
//     while (std::getline(file, line)) {
//         if (!line.empty() && line.back() == '\r') line.pop_back();
//         if (line.empty()) continue;
//         std::stringstream ss(line);
//         std::string cell;
//         std::vector<double> vals;
//         try {
//             while (std::getline(ss, cell, ',')) {
//                 // 去除可能存在的 Windows 回车符
//                 cell.erase(std::remove(cell.begin(), cell.end(), '\r'), cell.end()); 
//                 if (!cell.empty()) vals.push_back(std::stod(cell));
//             }
//         } catch (...) { continue; } 
        
//         // 兼容 6列(纯IMU) 或 7列(带时间) 格式，取前6列
//         if (vals.size() >= 6) {
//             raw_gyro.push_back(Vector3d(vals[0], vals[1], vals[2]));
//             raw_acc.push_back(Vector3d(vals[3], vals[4], vals[5]));
//         }
//     }

//     if (raw_gyro.empty()) return sim_data;

//     // 2. 加速度计归一化 (Auto-Scaling)
//     Vector3d acc_sum = Vector3d::Zero();
//     for (const auto& a : raw_acc) acc_sum += a;
//     Vector3d acc_mean_raw = acc_sum / raw_acc.size();
    
//     double norm_meas = acc_mean_raw.norm();
//     double scale_ratio = (norm_meas > 1e-5) ? target_g / norm_meas : 1.0; 
    
//     // 3. 构建全量 IMUData (用于 Nav)
//     std::vector<IMUData>& all_data = sim_data.data_nav;
//     double t_curr = 0.0;
//     for (size_t i = 0; i < raw_gyro.size(); ++i) {
//         IMUData d;
//         d.wm = raw_gyro[i] * ts;          
//         d.vm = (raw_acc[i] * scale_ratio) * ts; 
//         d.t = t_curr;
//         all_data.push_back(d);
//         t_curr += ts;
//     }
//     sim_data.time_nav = t_curr;

//     // 4. 数据切片 (Slicing) - 修正版
//     // 策略：粗对准用前 T1 秒，精对准用前 T2 秒 (包含 T1，KF通常需要从头收敛)，导航用 T2 之后的数据
    
//     size_t n_coarse = static_cast<size_t>(req_t_coarse / ts);
//     if (n_coarse > all_data.size()) n_coarse = all_data.size();
//     sim_data.data_coarse.assign(all_data.begin(), all_data.begin() + n_coarse);
//     sim_data.time_coarse = n_coarse * ts;

//     size_t n_fine = static_cast<size_t>(req_t_fine / ts); // 精对准通常包含粗对准段时间
//     if (n_fine > all_data.size()) n_fine = all_data.size();
//     sim_data.data_fine.assign(all_data.begin(), all_data.begin() + n_fine);
//     sim_data.time_fine = n_fine * ts;
    
//     // 导航数据：从精对准结束那一刻开始，直到文件结束
//     if (n_fine < all_data.size()) {
//         sim_data.data_nav.assign(all_data.begin() + n_fine, all_data.end());
//         sim_data.time_nav = sim_data.data_nav.size() * ts;
//     } else {
//         sim_data.data_nav.clear();
//         sim_data.time_nav = 0;
//     }

//     return sim_data;
// }