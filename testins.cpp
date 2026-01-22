#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <memory>
#include <Eigen/Dense>

#include "ins_wrapper.h"  

using namespace std;
using namespace Eigen;

std::vector<IMUData> load_and_process_data(const std::string& filename, double ts, double g_ref);

int main() {
    GLV glv; 
    string fog_path = "../fog.csv"; 
    double ts = 1.0 / 400.0;
    auto all_data = load_and_process_data(fog_path, ts, glv.g);
    if (all_data.empty()) return -1;

    // --- 配置对准参数 ---
    AlignConfig align_cfg;
    align_cfg.att_ref << -0.03258, 0.20927, 0.62977; 
    align_cfg.att_ref *= glv.deg;  
    align_cfg.pos_ref << 32.0286 * glv.deg, 118.8533 * glv.deg, 17.0; 
    align_cfg.phi_init_err << 0.01, 0.01, 0.1;
    align_cfg.phi_init_err *= glv.deg; 
    align_cfg.wvn_err << 0.001, 0.001, 0.001;
    align_cfg.eb_sigma = 0.2 * glv.dph;         
    align_cfg.db_sigma = 100 * glv.ug;          
    align_cfg.web_psd  = 0.021 * glv.dpsh;      
    align_cfg.wdb_psd  = 10 * glv.ugpsHz;

    // --- 全局时钟与控制 ---
    double global_T = 0.0;
    double align_duration = 300.0;

    // 实例化两个引擎
    AlignmentEngine aligner;
    aligner.Init(ts, glv, align_cfg);

    INSState ins; // 此时为空，稍后初始化
    Earth eth(glv);

    bool is_aligning = true;

    ofstream out_file("res_pure_global_timer.csv");
    out_file << "t,lat,lon,h,vE,vN,vU,roll,pitch,yaw,ebX,ebY,ebZ,dbX,dbY,dbZ" << endl;

    cout << "[System] Starting Simulation with Global Clock T..." << endl;

    // ============================================
    // 统一的主循环：每次循环代表一个 ts 的流逝
    // ============================================
    for (size_t k = 0; k < all_data.size(); ++k) {
        const auto& d = all_data[k];
        global_T = d.t; // 或者 global_T += ts;

        if (is_aligning) {
            // [阶段 1] 对准模式
            aligner.Step(d);

            // 检查是否完成
            if (global_T >= align_duration) {
                cout << "\n[System] Switching from Alignment to Navigation at T=" << global_T << endl;
                
                // 1. 获取对准结果
                AlignResult res = aligner.GetResult();
                
                // 2. 初始化惯导状态
                // 注意：这里用当前位置或参考位置均可，纯惯导通常用参考位置或对准结束位置
                ins = INSState(res.att, Vector3d::Zero(), align_cfg.pos_ref, ts, eth);
                ins.set_bias(res.eb, res.db);
                ins.is_align = false;
                ins.vertical_damping_mode = 1;
                
                // 3. 更新地球参数到当前位置
                eth.eupdate(ins.pos, ins.vn);

                // 4. 打印对准结果
                Vector3d eb_dph = res.eb * glv.rad * 3600.0;
                Vector3d db_ug  = res.db / glv.ug;
                cout << "  Est eb (deg/h): " << eb_dph.transpose() << endl;
                cout << "  Est db (ug):    " << db_ug.transpose() << endl;

                is_aligning = false; 
            }
        } 
        else {
            // [阶段 2] 纯惯导模式
            // INSState::update 已经是以 ts 为步长，直接调用
            ins.update(d.wm, d.vm, glv, eth);
            
            // 记录数据 (降频记录，例如每 0.1s 记一次)
            if (k % 40 == 0) {
                Vector3d att_deg = ins.att * glv.rad; 
                Vector3d pos_deg = ins.pos; 
                pos_deg(0) *= glv.rad; 
                pos_deg(1) *= glv.rad; 
                Vector3d eb_dph = ins.eb * glv.rad * 3600.0;
                Vector3d db_ug  = ins.db / glv.ug;

                out_file << fixed << setprecision(9) 
                         << global_T << ","
                         << pos_deg(0) << "," << pos_deg(1) << "," << ins.pos(2) << ","
                         << ins.vn(0) << "," << ins.vn(1) << "," << ins.vn(2) << ","
                         << att_deg(1) << "," << att_deg(0) << "," << att_deg(2) << "," 
                         << eb_dph(0) << "," << eb_dph(1) << "," << eb_dph(2) << ","
                         << db_ug(0) << "," << db_ug(1) << "," << db_ug(2) << endl;
            }
        }
    }

    cout << " Navigation Finished." << endl;
    return 0;
}



// [修改 1] 参数名统一为 g_ref
std::vector<IMUData> load_and_process_data(const std::string& filename, double ts, double g_ref) {
    std::vector<Vector3d> raw_gyro, raw_acc;
    std::ifstream file(filename);
    std::string line;
    
    if (!file.is_open()) { cerr << "Error: Cannot open " << filename << endl; return {}; }
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string cell;
        std::vector<double> vals;
        try {
            while (std::getline(ss, cell, ',')) {
                if (!cell.empty()) vals.push_back(std::stod(cell));
            }
        } catch (...) { continue; } 
        
        if (vals.size() >= 6) {
            Vector3d w(vals[0], vals[1], vals[2]);
            Vector3d a(vals[3], vals[4], vals[5]);
            raw_gyro.push_back(w);
            raw_acc.push_back(a);
        }
    }

    if (raw_gyro.empty()) return {};

    Vector3d acc_sum = Vector3d::Zero();
    for (const auto& a : raw_acc) acc_sum += a;
    Vector3d acc_mean_raw = acc_sum / raw_acc.size();
    
    double norm_meas = acc_mean_raw.norm();
    double scale_ratio = g_ref / norm_meas; 

    cout << "\n========== [Gravity Normalization] ==========" << endl;
    cout << "Ref Gravity (glv.g0):   " << fixed << setprecision(8) << g_ref << " m/s^2" << endl;
    cout << "Meas Acc Norm (Raw):    " << norm_meas << " m/s^2" << endl;
    cout << "Scale Ratio:            " << scale_ratio << endl;
    cout << "=============================================" << endl;

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
