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
#include "cai_sim.h"

using namespace std;
using namespace Eigen;
using namespace cai;

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


class IMUChainedLoader {
public:
    // 构造函数：传入所有要读的文件名列表
    IMUChainedLoader(const std::vector<std::string>& file_list, double ts, double target_g) 
        : files_(file_list), ts_(ts), target_g_(target_g), 
          current_file_idx_(0), current_data_idx_(0) {
        
        // 初始化时预加载第一个文件
        if (!files_.empty()) {
            LoadNextFile();
        }
    }

    // 核心接口：获取下一个数据点
    // 返回 true 表示成功拿到数据
    // 返回 false 表示所有文件都读完了
    bool Next(IMUData& out_data) {
        // 1. 如果当前 buffer 读完了，尝试加载下一个文件
        if (current_data_idx_ >= current_buffer_.size()) {
            // 如果还有文件没读，加载下一个
            if (current_file_idx_ < files_.size()) {
                LoadNextFile();
                // 如果加载完发现还是空的（可能是空文件），递归尝试再下一个，或者返回 false
                if (current_buffer_.empty()) return Next(out_data); 
            } else {
                // 所有文件都读完了
                return false; 
            }
        }

        // 2. 取出数据
        out_data = current_buffer_[current_data_idx_];
        
        // 3. 指针后移
        current_data_idx_++;
        return true;
    }

    // 获取总体进度（可选，用于显示）
    std::string GetCurrentFileName() const {
        if (current_file_idx_ > 0 && current_file_idx_ <= files_.size()) {
            return files_[current_file_idx_ - 1];
        }
        return "None";
    }

private:
    // 内部函数：释放旧内存，加载新文件
    void LoadNextFile() {
        // 1. 彻底释放旧内存 (你的核心诉求)
        std::vector<IMUData>().swap(current_buffer_); 
        current_data_idx_ = 0;

        // 2. 获取文件名
        std::string filename = files_[current_file_idx_];
        std::cout << "[Loader] Switching to file: " << filename << " ..." << std::endl;

        // 3. 调用你提供的 LoadIMUData 函数
        current_buffer_ = LoadIMUData(filename, ts_, target_g_);
        
        std::cout << "[Loader] Loaded " << current_buffer_.size() << " epochs. Memory Usage Refreshed." << std::endl;

        // 4. 文件索引 +1
        current_file_idx_++;
    }

    // 成员变量
    std::vector<std::string> files_;
    std::vector<IMUData> current_buffer_; // 当前驻留内存的数据块
    size_t current_file_idx_; // 下一个要读的文件索引
    size_t current_data_idx_; // 当前 buffer 读到哪里了
    
    double ts_;
    double target_g_;
};


Vector3d LocalToGeo(const Vector3d& xyz_enu, const Vector3d& pos_ref_rad, const GLV& glv) {
    double lat = pos_ref_rad(0); double h = pos_ref_rad(2);
    double sl = sin(lat); double cl = cos(lat); 
    double e2 = glv.e2; double Re = glv.Re;
    double sq = sqrt(1.0 - e2 * sl * sl);
    double RN = Re / sq; double RM = RN * (1.0 - e2) / (1.0 - e2 * sl * sl); 
    return Vector3d(lat + xyz_enu(1)/(RM+h), pos_ref_rad(1) + xyz_enu(0)/((RN+h)*cl), h + xyz_enu(2));
}


struct AtomicAlignResult {
    Matrix3d Cnb;
    Vector3d att; // rad
    Vector3d db;  // 加计零偏估计 (m/s^2)
};

AtomicAlignResult Align_Atomic_6H(const std::vector<IMUData>& fog_data, 
                                  const AtomicGyroSimulator& atom_sim,
                                  const Earth& eth,
                                  double ts) {
    cout << "\n[AtomicAlign] Generating 6h atomic data & Running independent alignment..." << endl;
    
    // 1. 累积器
    Vector3d sum_wb_atom = Vector3d::Zero(); // 原子陀螺角速度累积
    Vector3d sum_fb_real = Vector3d::Zero(); // 真实加计比力累积
    
    // 仿真参数
    double T_cycle = 2.0; // 原子更新周期 2s
    int samples_per_cycle = static_cast<int>(T_cycle / ts); // 800
    
    // 临时累积器 (用于把加计降采样到 2s)
    Vector3d chunk_fb_sum = Vector3d::Zero();
    int chunk_cnt = 0;
    int atom_epochs = 0;

    // 复制一个仿真器实例，避免影响主程序的随机数序列状态（虽然影响也不大）
    AtomicGyroSimulator atom_runner = atom_sim; 

    // 2. 遍历 6小时数据
    for (const auto& epoch : fog_data) {
        chunk_fb_sum += epoch.vm / ts; // m/s^2
        chunk_cnt++;
        
        if (chunk_cnt >= samples_per_cycle) {
            // A. 生成这一帧的原子数据 (模拟原子陀螺输出)
            Vector3d wb_atom = atom_runner.Measure(); 
            
            // B. 获取这一帧的真实加计数据 (2s 平均值)
            Vector3d fb_real = chunk_fb_sum / samples_per_cycle;
            
            // C. 全局累积
            sum_wb_atom += wb_atom;
            sum_fb_real += fb_real;
            atom_epochs++;
            
            // 重置分块
            chunk_fb_sum.setZero();
            chunk_cnt = 0;
        }
    }
    
    if (atom_epochs == 0) return {};

    // 3. 计算 6小时 全局均值
    Vector3d wb_mean = sum_wb_atom / atom_epochs; // b系 原子角速度均值
    Vector3d fb_mean = sum_fb_real / atom_epochs; // b系 真实比力均值
    
    cout << "  Processed " << atom_epochs << " atomic epochs." << endl;
    
    // 4. 双矢量定姿 (Triad Algorithm)
    //    n 系基准: 重力方向(-gn) 和 地球自转(wien)
    Vector3d vn1 = -eth.gn; 
    Vector3d vn2 = eth.wien;
    
    //    b 系测量: 加计均值 和 原子均值
    Vector3d vb1 = fb_mean;
    Vector3d vb2 = wb_mean;
    
    //    构造旋转矩阵 Cnb
    //    主矢量选用重力(加计)，因为加计信噪比通常比陀螺高，且此时原子陀螺也是基于重力对齐的姿态生成的
    Vector3d n_e1 = vn1.normalized();
    Vector3d n_e2 = vn1.cross(vn2).normalized(); // East
    Vector3d n_e3 = n_e1.cross(n_e2).normalized(); // North
    Matrix3d Mn; Mn << n_e1, n_e2, n_e3;
    
    Vector3d b_e1 = vb1.normalized();
    Vector3d b_e2 = vb1.cross(vb2).normalized();
    Vector3d b_e3 = b_e1.cross(b_e2).normalized();
    Matrix3d Mb; Mb << b_e1, b_e2, b_e3;
    
    Matrix3d Cnb = Mn * Mb.transpose();
    Vector3d att = INSMath::m2att(Cnb);
    
    return {Cnb, att};
}

// ==========================================
// [工具] 纯原子导航验证函数 (完整实现)
// ==========================================
void Run_Pure_Atomic_Nav(const Vector3d& pos_ref, const Vector3d& att_initial, 
                         AtomicGyroSimulator& atom_gyro, 
                         AtomicAccSimulator& atom_acc,
                         double ts) {
    cout << "\n==================================================" << endl;
    cout << "           PURE ATOMIC NAVIGATION TEST            " << endl;
    cout << "==================================================" << endl;
    
    // 1. 初始化引擎
    SinsEngine test_eng(ts);
    // 必须手动初始化 Earth 和 INS
    test_eng.eth.update(pos_ref, Vector3d::Zero()); 
    test_eng.ins = INSState(att_initial, Vector3d::Zero(), pos_ref, ts, test_eng.eth); 
    
    cout << " [Init] Position: " << pos_ref.transpose() << endl;
    cout << " [Init] Attitude: " << att_initial.transpose() << endl;

    double total_time = 24.0 * 3600.0;
    size_t steps = static_cast<size_t>(total_time / ts);
    double max_drift = 0;
    GLV glv;

    // 2. 跑循环
    for (size_t i = 0; i < steps; ++i) {
        // A. 生成数据 
        Vector3d wm = atom_gyro.Measure() * ts; 
        Vector3d vm = atom_acc.Measure() * ts;  

        // B. 惯导解算 (纯机械编排)
        test_eng.Step_Nav(wm, vm);
        
        // C. 记录漂移 (每1小时打印一次)
        if (i > 0 && i % (3600 * 400) == 0) { 
            double lat_err = (test_eng.ins.pos(0) - pos_ref(0)) * glv.Re;
            double lon_err = (test_eng.ins.pos(1) - pos_ref(1)) * glv.Re * cos(pos_ref(0));
            double drift = sqrt(lat_err*lat_err + lon_err*lon_err);
            
            if (drift > max_drift) max_drift = drift;
            
            cout << fixed << setprecision(3) 
                 << " Time: " << (i*ts)/3600.0 << "h | Drift: " << drift << " m" << endl;
        }
    }
    cout << "--------------------------------------------------" << endl;
    cout << " [Result] Pure Atomic Max Drift (24h): " << max_drift << " m" << endl;
    cout << "==================================================\n" << endl;
}



#endif // SUPPORT_H

