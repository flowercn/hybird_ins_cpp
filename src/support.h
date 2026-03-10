#ifndef SUPPORT_H
#define SUPPORT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <random>
#include <Eigen/Dense>
#include "../psins/sins_engine.h"  
#include "cai_sim.h"

using namespace std;
using namespace Eigen;
using namespace cai;

// ==========================================
// 用于 FGO 解析回溯的日志缓存结构
// ==========================================
struct LogPoint {
    double nav_time;
    Vector3d att;
    Vector3d vn;
    Vector3d pos;
    Matrix3d Cnb;
    MatrixXd J_at_t;
    
    LogPoint(double t, const Vector3d& a, const Vector3d& v, const Vector3d& p, const Matrix3d& C, const MatrixXd& j) 
        : nav_time(t), att(a), vn(v), pos(p), Cnb(C), J_at_t(j) {}
};

enum class AlgoType {
    PURE_FOG,   // 纯光纤机械编排 (用于提取最大发散误差)
    ES_FGO,     // 因子图解析平滑 (O(1) 线性映射)
    EKF,        // 传统重积分滤波
    TRUE_FGO    // 真正的非线性迭代因子图优化 (Gauss-Newton + 预积分)
};

// 实验配置参数
struct ExperimentConfig {
    string name;           // 实验名
    AlgoType algo;         // 算法类型
    double duration_hours; // 仿真运行时间 (小时)
    double t_active;       // 有效时间 (s)
    double t_dead;         // 死区时间 (s)
    bool use_atomic_acc;   // 是否使用原子加计
    string output_file;    // 输出文件名
    
    // ====== 【新增】动态时变漂移参数 ======
    bool inject_time_varying_drift = false; // 是否注入动态漂移
    double drift_amplitude = 0.0;           // 漂移幅值 (°/h)
    double drift_period = 1.0;              // 漂移周期 (h)
};

// 共享的上下文环境 (避免重复传参)
struct SimulationContext {
    double ts;
    GLV glv;
    Vector3d pos_ref;
    Earth eth;
    double local_g;
    std::vector<std::string> file_list; // 数据文件列表
    
    // 核心基准：对准后的姿态和零偏
    Vector3d att_align; 
    Vector3d eb_align;
    Vector3d db_align;
    Matrix3d Cnb_align; // 【关键】锁定的姿态基准
};

// ==============================================================================
// 辅助函数：日志记录
// ==============================================================================
inline void LogData(ofstream& log, double time, const INSState& ins, const SimulationContext& ctx, double drift, double true_eb_x_dph = 0.0) {
    log << fixed << setprecision(9) << time << ","
        << (ins.pos(0) - ctx.pos_ref(0)) * ctx.glv.Re << ","
        << (ins.pos(1) - ctx.pos_ref(1)) * ctx.glv.Re * cos(ctx.pos_ref(0)) << ","
        << drift << ","
        << ins.vn(0) << ","
        << ins.vn(1) << ","
        << ins.vn(2) << ","
        << ins.att(0) / ctx.glv.deg << ","
        << ins.att(1) / ctx.glv.deg << ","
        << ins.att(2) / ctx.glv.deg << ","
        << ins.eb(0) / ctx.glv.deg * 3600.0 << ","
        << ins.eb(1) / ctx.glv.deg * 3600.0 << ","
        << ins.eb(2) / ctx.glv.deg * 3600.0 << ","
        << ins.db(0) / ctx.glv.ug << ","
        << ins.db(1) / ctx.glv.ug << ","
        << ins.db(2) / ctx.glv.ug << ","
        << true_eb_x_dph << "\n";
}

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

#endif // SUPPORT_H
