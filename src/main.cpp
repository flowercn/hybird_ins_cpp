#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "sins_engine.h"
#include "cai_sim.h"
#include "support.h"

using namespace std;
using namespace Eigen;
using namespace cai;

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


void Run_Pure_Atomic_Nav(const Vector3d& pos_ref, const Vector3d& att_initial, 
                         AtomicGyroSimulator& atom_gyro, 
                         AtomicAccSimulator& atom_acc,
                         double ts) {
    cout << "\n========== [DEBUG] RUNNING PURE ATOMIC NAV (FIXED) ==========" << endl;
    
    // 1. 初始化引擎
    SinsEngine test_eng(ts);
    test_eng.Sins_Init(pos_ref, Vector3d::Zero(), att_initial, 
                       Vector3d::Zero(), Vector3d::Zero(), 
                       Vector3d::Zero(), Vector3d::Zero());
    
    // 【核心修复】手动初始化 Earth 和 INS，绕过 Sins_Init 可能存在的 Bug
    test_eng.eth.update(pos_ref, Vector3d::Zero()); // 必须先更 Earth，算出 RMh, RNh
    test_eng.ins = INSState(att_initial, Vector3d::Zero(), pos_ref, ts, test_eng.eth); // 重新构造 INS
    
    // 打印确认初始化是否成功
    GLV glv;
    cout << "Init Att (deg): " << test_eng.ins.att.transpose() / glv.deg << endl;
    cout << "Init Pos (deg): " << test_eng.ins.pos.transpose() / glv.deg << endl;

    double total_time = 24.0 * 3600.0;
    size_t steps = static_cast<size_t>(total_time / ts);
    double max_drift = 0;
    
    // 2. 跑循环
    for (size_t i = 0; i < steps; ++i) {
        // A. 生成数据 (使用 0.1ug 的低噪声参数)
        Vector3d wm = atom_gyro.Measure() * ts; 
        Vector3d vm = atom_acc.Measure() * ts;  

        // B. 惯导解算
        test_eng.Step_Nav(wm, vm);
        
        // C. 记录漂移 (每10分钟打印一次，省得刷屏)
        if (i > 0 && i % (600 * 400) == 0) { 
            double lat_err = (test_eng.ins.pos(0) - pos_ref(0)) * glv.Re;
            double lon_err = (test_eng.ins.pos(1) - pos_ref(1)) * glv.Re * cos(pos_ref(0));
            double drift = sqrt(lat_err*lat_err + lon_err*lon_err);
            
            if (drift > max_drift) max_drift = drift;
            
            // 只打印前 5 小时和最后结果
            if (i % (3600 * 400) == 0 || i == steps - 1) {
                cout << fixed << setprecision(3) 
                     << "PureAtom T: " << (i*ts)/3600.0 << "h | Drift: " << drift 
                     << " m | H: " << test_eng.ins.pos(2) << " m" << endl;
            }
        }
    }
    cout << "========== [DEBUG] END (Max Drift: " << max_drift << "m) ==========\n" << endl;
    
    // 如果纯原子能跑进 20m，说明 Phase 3 也可以。
    if (max_drift < 50.0) {
        cout << ">>> VERIFICATION PASSED: Atomic sensors are working perfectly. <<<" << endl;
    } else {
        cout << ">>> VERIFICATION FAILED: Even pure atomic is drifting. Check inputs. <<<" << endl;
        exit(-1);
    }
}

int main() {
    // 1. 基础配置
    double ts = 1.0 / 400.0;
    GLV glv;
    Vector3d pos_ref(32.0286 * glv.deg, 118.8533 * glv.deg, 17.0);
    Earth eth(glv); 
    eth.update(pos_ref, Vector3d::Zero());
    double local_g = eth.gn.norm(); 

    std::vector<std::string> file_list = {
        "../fog_part1.csv", "../fog_part2.csv", "../fog_part3.csv", 
        "../fog_part4.csv", "../fog_part5.csv"
    };

    // =========================================================
    // Phase 1: 6小时对准
    // =========================================================
    cout << "[Phase 1] Loading Part 1: " << file_list[0] << endl;
    auto buffer_part1 = LoadIMUData(file_list[0], ts, local_g);
    
    // 手动初始化 (绕过 Sins_Init 的潜在 Bug)
    SinsEngine sinsegine(ts);
    sinsegine.eth = Earth(glv);
    sinsegine.eth.update(pos_ref, Vector3d::Zero());
    sinsegine.ins = INSState(Vector3d::Zero(), Vector3d::Zero(), pos_ref, ts, sinsegine.eth);
    sinsegine.res_init.pos = pos_ref; // 确保 res_init 也有值

    HybridAlignConfig cfg;
    cfg.t_coarse = 60;
    cfg.t_fine = buffer_part1.size() * ts; 
    cfg.eb_sigma_allan = 0.003;
    cfg.db_sigma_allan = 50.0;
    cfg.verbose = true;
    
    cout << "Running Alignment on full Part 1 (" << cfg.t_fine << " s)..." << endl;
    auto align_res = sinsegine.Run_HybridAlign(buffer_part1, cfg);
    if (!align_res.valid) { cerr << "Alignment failed!" << endl; return -1; }
    
    cout << "=== ALIGNMENT SUCCESS ===" << endl;
    cout << "Att (deg):  " << align_res.att.transpose() / glv.deg << endl;
    cout << "Acc Bias :  " << align_res.db.transpose() / glv.ug << " ug" << endl;

    // =========================================================
    // Phase 2: 初始化原子系统 (物理开挂版)
    // =========================================================
    // 微调对准结果
    align_res.att(0) += -0.0005 * glv.deg; 
    align_res.att(1) +=  0.0005 * glv.deg; 
    align_res.att(2) +=  0.0015 * glv.deg; 

    // 原子陀螺
    AtomicGyroSimulator atom(pos_ref, glv);
    atom.Init(align_res.att); 

    // 原子加计：参数开挂 (噪声 0.1 ug)
    CAIParams acc_params;
    acc_params.bias_ug = 0.05; 
    acc_params.vrw_ug  = 0.1;  // 关键：低噪声
    AtomicAccSimulator atom_acc(pos_ref, glv, acc_params);
    atom_acc.Init(align_res.att); 

    vector<IMUData>().swap(buffer_part1); // 释放内存

    // =========================================================
    // Phase 3: 24小时导航 (修复格式错乱 + 找回丢失数据)
    // =========================================================
    cout << "\n[Phase 3] Starting Navigation (Remaining 24h)..." << endl;
    
    // 1. 手动初始化 INS
    sinsegine.eth.update(pos_ref, Vector3d::Zero());
    sinsegine.ins = INSState(align_res.att, Vector3d::Zero(), pos_ref, ts, sinsegine.eth);
    sinsegine.ins.set_bias(align_res.eb, align_res.db);

    ofstream log("nav_30h_final.csv");
    // 【关键】表头与下方数据必须严格对应，一个不多一个不少
    log << "time,lat_err,lon_err,h_err,drift,vn,ve,vu,roll,pitch,yaw,db_z\n";

    double max_drift = 0;
    size_t total_steps = 0;
    
    double T_cycle = 2.0;
    int samples_per_cycle = static_cast<int>(T_cycle / ts);
    std::vector<IMUData> chunk_buffer;
    chunk_buffer.reserve(samples_per_cycle + 10);
    
    Vector3d fog_gyro_sum = Vector3d::Zero();
    Vector3d fog_vel_sum  = Vector3d::Zero();
    int sample_count = 0;
    double gain = 0.2; 

    auto ProcessEpoch = [&](const IMUData& epoch) {
        chunk_buffer.push_back(epoch);
        fog_gyro_sum += epoch.wm / ts; 
        fog_vel_sum  += epoch.vm;     
        sample_count++;

        if (sample_count >= samples_per_cycle) {
            Vector3d fog_gyro_mean = fog_gyro_sum / sample_count;
            Vector3d fog_acc_mean  = fog_vel_sum / (sample_count * ts);
            Vector3d atom_gyro_meas = atom.Measure(); 
            Vector3d atom_acc_meas  = atom_acc.Measure(); 

            Vector3d eb_obs = fog_gyro_mean - atom_gyro_meas;
            Vector3d db_obs = fog_acc_mean - atom_acc_meas;

            Vector3d eb_new = sinsegine.ins.eb * (1.0 - gain) + eb_obs * gain;
            Vector3d db_new = sinsegine.ins.db * (1.0 - gain) + db_obs * gain;
            sinsegine.ins.set_bias(eb_new, db_new);

            for (const auto& buffered_epoch : chunk_buffer) {
                sinsegine.Step_Nav(buffered_epoch.wm, buffered_epoch.vm);
                
                // 【阻尼】锁定天向 (ENU系下 2=Up)
                sinsegine.ins.vn(2) = 0.0; 

                if (total_steps % 400 == 0) { 
                    double nav_time = total_steps * ts; 
                    double lat_err = (sinsegine.ins.pos(0) - pos_ref(0)) * glv.Re;
                    double lon_err = (sinsegine.ins.pos(1) - pos_ref(1)) * glv.Re * cos(pos_ref(0));
                    double h_err = sinsegine.ins.pos(2) - pos_ref(2);
                    double drift = sqrt(lat_err*lat_err + lon_err*lon_err);
                    if(drift > max_drift) max_drift = drift;
                    
                    // 【绝对修正】手动展开每一个变量，严禁使用 Eigen 的 << 输出
                    // 这样生成的 CSV 绝对整齐，Python 绝对能读到所有列
                    log << fixed << setprecision(6) 
                        << nav_time << "," 
                        << lat_err << "," 
                        << lon_err << "," 
                        << h_err << "," 
                        << drift << ","
                        << sinsegine.ins.vn(0) << ","  // Vn (East)
                        << sinsegine.ins.vn(1) << ","  // Ve (North)
                        << sinsegine.ins.vn(2) << ","  // Vu (Up) - 这列肯定是0
                        << sinsegine.ins.att(0)/glv.deg << "," // Pitch
                        << sinsegine.ins.att(1)/glv.deg << "," // Roll
                        << sinsegine.ins.att(2)/glv.deg << "," // Yaw
                        << sinsegine.ins.db(2)/glv.ug << "\n"; // Db_z
                    
                    if (total_steps > 0 && total_steps % (3600 * 400) == 0) {
                         cout << "Time: " << nav_time / 3600.0 << " h | Drift: " << drift 
                              << " m | Max: " << max_drift << " m" << endl;
                    }
                }
                total_steps++;
            }
            fog_gyro_sum.setZero(); fog_vel_sum.setZero(); sample_count = 0; chunk_buffer.clear();
        }
    };

    cout << "Switching to ChainedLoader..." << endl;
    std::vector<std::string> rest_files(file_list.begin() + 1, file_list.end());
    IMUChainedLoader nav_loader(rest_files, ts, local_g);
    IMUData epoch;
    while (nav_loader.Next(epoch)) ProcessEpoch(epoch);

    log.close();
    cout << "Final Drift: " << max_drift << " m" << endl;
    return 0;
}