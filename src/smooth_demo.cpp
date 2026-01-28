#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <random> 
#include "sins_engine.h"
#include "support.h"

using namespace std;
using namespace Eigen;

// =================================================================================
// 物理仿真核心：原子陀螺(CAIG) 模拟器
// 物理模型：原子陀螺测量的是"平台真实运动" = -(INS计算的姿态变化)
// 在稳定平台系统中，INS姿态变化被反向控制到平台上
// =================================================================================
class AtomicGyroSimulator {
public:
    double arw_psd;        // 角度随机游走 PSD (rad/sqrt(s))
    double ts;             // 单步时间
    mt19937 rng;
    Vector3d noise_accum;  // 累积的ARW噪声

    AtomicGyroSimulator(double arw_rad_sqrt_s, double sample_ts, const GLV& glv) 
        : arw_psd(arw_rad_sqrt_s), ts(sample_ts), rng(42) {
        noise_accum.setZero();
        cout << "[AtomGyro] Initialized. ARW = " << arw_rad_sqrt_s / glv.deg * 60.0 
             << " deg/sqrt(h)" << endl;
    }

    // 每步积累ARW噪声
    void Step() {
        double sigma = arw_psd * sqrt(ts);
        normal_distribution<double> dist(0.0, sigma);
        noise_accum += Vector3d(dist(rng), dist(rng), dist(rng));
    }
    
    // 给定INS某时刻的姿态，返回原子陀螺"测量"的平台姿态
    // 平台姿态 = INS姿态的逆（因为反向控制）+ 累积噪声
    Quaterniond GetPlatformAttitude(const Quaterniond& q_ins) const {
        Quaterniond q_noise = INSMath::rv2q(noise_accum);
        // 平台姿态 = q_ins^{-1}，然后叠加噪声
        return q_ins.inverse() * q_noise;
    }
};

// =================================================================================
// 核心算法：回溯平滑 (Backtracking Smoothing)
// 
// 物理原理：
// - INS计算姿态变化 dq_ins = q_161^{-1} * q_800
// - 原子测量平台变化 dq_atom = (q_161^{-1})^{-1} * (q_800^{-1}) = q_161 * q_800^{-1} = dq_ins^{-1}
// - 所以理想情况：dq_ins * dq_atom = I
// - 如果有误差：dq_ins * dq_atom = dq_error ≠ I
// - 这个误差反映了零偏估计不准导致的姿态累积误差
// =================================================================================
void Run_Backtracking_Smoothing(const SinsEngine& engine_template,
                                const std::vector<IMUData>& data_segment,
                                INSState& state_curr, 
                                AtomicGyroSimulator& atom_sim,
                                std::vector<INSState>& segment_traj) 
{
    if (data_segment.size() < 800) {
        cerr << "[Error] Segment size < 800!" << endl;
        return;
    }
    
    const int IDX_161 = 160;  // 第161个点，索引160
    const int IDX_800 = 799;  // 第800个点，索引799
    
    // 1. 试跑 INS 800步，同时推进原子陀螺噪声累积
    SinsEngine trial_engine = engine_template; 
    trial_engine.ins = state_curr;
    
    Quaterniond q_ins_161, q_ins_800;

    for (int i = 0; i < 800; ++i) {
        trial_engine.Step_Nav(data_segment[i].wm, data_segment[i].vm);
        atom_sim.Step();  // 累积ARW噪声
        
        if (i == IDX_161) q_ins_161 = trial_engine.GetQnb();
        if (i == IDX_800) q_ins_800 = trial_engine.GetQnb();
    }

    // 2. 获取原子陀螺"测量"的平台姿态
    Quaterniond q_atom_161 = atom_sim.GetPlatformAttitude(q_ins_161);
    Quaterniond q_atom_800 = atom_sim.GetPlatformAttitude(q_ins_800);

    // 3. 计算相对姿态增量
    // INS: dq_ins = q_161^{-1} * q_800  （载体从161转到800的增量）
    // Atom: dq_atom = q_atom_161^{-1} * q_atom_800 （平台从161转到800的增量）
    Quaterniond dq_ins = q_ins_161.inverse() * q_ins_800;
    Quaterniond dq_atom = q_atom_161.inverse() * q_atom_800;
    
    // 4. 计算误差
    // 理想情况：dq_atom = dq_ins^{-1}，即 dq_ins * dq_atom = I
    // 实际：dq_ins * dq_atom = dq_error
    Quaterniond dq_error = dq_ins * dq_atom;
    Vector3d phi_error = INSMath::q2rv(dq_error);

    // 5. 计算零偏修正量
    double dt_interval = (IDX_800 - IDX_161) * engine_template.ins.ts;
    
    // 误差分析：
    // 如果 eb_est 偏小，INS多转了，平台也多转了，原子测到更大的负角度
    // dq_ins > I, dq_atom 更负，dq_ins * dq_atom 的符号...
    // 需要增大 eb_est
    Vector3d delta_bg = phi_error / dt_interval; 
    state_curr.eb += delta_bg;

    // 6. 用修正后的零偏重积分
    segment_traj.clear();
    segment_traj.reserve(800);
    
    SinsEngine final_engine = engine_template;
    final_engine.ins = state_curr;
    
    for (int i = 0; i < 800; ++i) {
        final_engine.Step_Nav(data_segment[i].wm, data_segment[i].vm);
        segment_traj.push_back(final_engine.ins);
    }
    
    state_curr = final_engine.ins;
}

// =================================================================================
// 主函数
// =================================================================================
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
    Vector3d att_truth(0.03404, 0.31606, 0.61788); // 你的真值
    engine.res_coarse.valid = true;
    engine.res_coarse.att = att_truth * glv.deg; // 注意单位
    engine.res_coarse.vel.setZero();
    engine.res_coarse.pos = pos_nj;
    engine.res_coarse.eb.setZero();
    engine.res_coarse.db.setZero();
    if (engine.res_coarse.valid) {
        Vector3d att = engine.res_coarse.att /glv.deg;
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
    // 3. 平滑导航（完整循环）
    // ---------------------------------------------------------
    cout << ">>> [3/3] Running Smoothed Navigation..." << endl;
    
    // DEBUG: 比较精对准后的零偏
    cout << "ins.eb (deg/h): " << (engine.ins.eb / glv.deg * 3600.0).transpose() << endl;
    cout << "res_fine.eb (deg/h): " << (engine.res_fine.eb / glv.deg * 3600.0).transpose() << endl;
    
    // [TEST] 不设置零偏，用 eb=0 作为基线
    // engine.ins.set_bias(engine.res_fine.eb, engine.res_fine.db);
    cout << "Using eb=0 for navigation baseline test" << endl;
    
    // 原子陀螺参数
    double arw_deg_sqrt_h = 2e-4;
    double arw_rad_sqrt_s = arw_deg_sqrt_h * glv.deg / 60.0;
    
    // 关键：继承精对准后的 engine.ins（包含正确的 wm_last, vm_last）
    // 测试时不设置零偏
    // engine.ins.set_bias(engine.res_fine.eb, engine.res_fine.db);
    
    INSState current_state = engine.ins;
    
    Vector3d p0 = engine.GetPosDeg();
    
    cout << "Initial pos (deg): " << p0.transpose() << endl;
    cout << "Initial att (deg): " << engine.GetAttDeg().transpose() << endl;
    cout << "Initial eb (deg/h): " << (current_state.eb / glv.deg * 3600.0).transpose() << endl;
    
    // 参数设置
    const int SEGMENT_SIZE = 800;
    const int IDX_161 = 160;
    const int IDX_800 = 799;
    double dt_161_800 = (IDX_800 - IDX_161) * ts;
    
    size_t nav_start_idx = 0;  // data_nav 已经是从精对准结束后的数据
    size_t total_nav_steps = sim.data_nav.size();
    size_t num_segments = total_nav_steps / SEGMENT_SIZE;
    
    cout << "Nav data points: " << total_nav_steps << endl;
    cout << "Number of 2s segments: " << num_segments << endl;
    
    // 输出文件
    ofstream f_bias("smooth_bias.csv");
    f_bias << "Time,Eb_X,Eb_Y,Eb_Z" << endl;
    ofstream f_traj("smooth_traj.csv");
    f_traj << "Time,Lat,Lon,H,Roll,Pitch,Yaw" << endl;
    
    // 随机数生成器（原子陀螺噪声）
    mt19937 rng(42);
    double sigma_total = arw_rad_sqrt_s * sqrt(dt_161_800);
    normal_distribution<double> dist(0.0, sigma_total);
    
    vector<INSState> full_traj;
    full_traj.reserve(num_segments * SEGMENT_SIZE);
    
    // ===========================================
    // 单步验证测试：检查之前99%结果的来源
    // ===========================================
    cout << "\n=== Single Segment Verification ===" << endl;
    {
        vector<IMUData> seg0;
        for (int i = 0; i < SEGMENT_SIZE; ++i) {
            seg0.push_back(sim.data_nav[i]);
        }
        
        INSState state0 = engine.ins;
        
        // === 第一次试跑 ===
        SinsEngine trial1 = engine;
        trial1.ins = state0;
        
        Quaterniond q1_161, q1_800;
        for (int i = 0; i < SEGMENT_SIZE; ++i) {
            trial1.Step_Nav(seg0[i].wm, seg0[i].vm);
            if (i == IDX_161) q1_161 = trial1.GetQnb();
            if (i == IDX_800) q1_800 = trial1.GetQnb();
        }
        
        Quaterniond dq_ins1 = q1_161.inverse() * q1_800;
        Vector3d phi_noise(dist(rng), dist(rng), dist(rng));
        Quaterniond dq_atom = dq_ins1.inverse() * INSMath::rv2q(phi_noise);
        
        // 残差定义1：dq_ins * dq_atom（这是噪声）
        Quaterniond dq_err1 = dq_ins1 * dq_atom;
        Vector3d phi_err1 = INSMath::q2rv(dq_err1);
        Vector3d delta_eb = phi_err1 / dt_161_800;
        
        // === 之前可能的错误逻辑：重跑后用【同一个 dq_atom】计算残差 ===
        INSState state_corrected = state0;
        state_corrected.eb = state0.eb + delta_eb;
        
        SinsEngine trial2 = engine;
        trial2.ins = state_corrected;
        
        Quaterniond q2_161, q2_800;
        for (int i = 0; i < SEGMENT_SIZE; ++i) {
            trial2.Step_Nav(seg0[i].wm, seg0[i].vm);
            if (i == IDX_161) q2_161 = trial2.GetQnb();
            if (i == IDX_800) q2_800 = trial2.GetQnb();
        }
        
        Quaterniond dq_ins2 = q2_161.inverse() * q2_800;
        
        // 错误方法：用同一个 dq_atom（没有更新）
        Quaterniond dq_err2_wrong = dq_ins2 * dq_atom;  // <-- 这里用的是旧的 dq_atom！
        Vector3d phi_err2_wrong = INSMath::q2rv(dq_err2_wrong);
        
        // 正确方法：重新仿真原子陀螺
        Quaterniond dq_atom2 = dq_ins2.inverse() * INSMath::rv2q(phi_noise);
        Quaterniond dq_err2_correct = dq_ins2 * dq_atom2;
        Vector3d phi_err2_correct = INSMath::q2rv(dq_err2_correct);
        
        cout << "phi_noise norm: " << phi_noise.norm() << " rad" << endl;
        cout << "phi_err1 (before): " << phi_err1.norm() << " rad" << endl;
        cout << "phi_err2 WRONG (same dq_atom): " << phi_err2_wrong.norm() << " rad" << endl;
        cout << "phi_err2 CORRECT (updated dq_atom): " << phi_err2_correct.norm() << " rad" << endl;
        cout << "\nRatio WRONG: " << phi_err2_wrong.norm() / phi_err1.norm() << endl;
        cout << "Ratio CORRECT: " << phi_err2_correct.norm() / phi_err1.norm() << endl;
        
        cout << "\n*** If ratio_wrong is ~1e-5, that explains the 99% reduction! ***" << endl;
        cout << "*** The bug was: dq_atom should change when INS changes in a stabilized platform! ***" << endl;
    }
    cout << "=== End Verification ===\n" << endl;
    
    // 零偏估计状态（使用简单的滑动平均）
    Vector3d eb_estimate = engine.res_fine.eb;
    const double alpha = 0.005;  // 滤波系数，值越小越平滑
    
    // DEBUG: 打印第一个和最后一个 IMU 数据
    if (!sim.data_nav.empty()) {
        cout << "First IMU data time: " << sim.data_nav[0].t << endl;
        cout << "Last IMU data time: " << sim.data_nav.back().t << endl;
        cout << "Total data_nav size: " << sim.data_nav.size() << endl;
    }
    
    // 简单版本：直接跑所有数据（不分段）
    for (size_t i = 0; i < sim.data_nav.size(); ++i) {
        engine.Step_Nav(sim.data_nav[i].wm, sim.data_nav[i].vm);
    }
    current_state = engine.ins;
    
    /*
    for (size_t seg = 0; seg < num_segments; ++seg) {
        size_t base_idx = nav_start_idx + seg * SEGMENT_SIZE;
        
        // 直接用 engine 跑所有数据（测试）
        for (int i = 0; i < SEGMENT_SIZE; ++i) {
            engine.Step_Nav(sim.data_nav[base_idx + i].wm, sim.data_nav[base_idx + i].vm);
        }
        
        // C. 更新 current_state
        current_state = engine.ins;
        
        // D. 记录零偏
        Vector3d eb_dph = current_state.eb / glv.deg * 3600.0;
        f_bias << sim.data_nav[base_idx].t << "," 
               << eb_dph.x() << "," << eb_dph.y() << "," << eb_dph.z() << endl;
        
        // K. 进度显示
        if ((seg + 1) % 100 == 0) {
            double t_now = sim.data_nav[base_idx + SEGMENT_SIZE - 1].t;
            cout << "    Processed " << (seg + 1) << " segments (" 
                 << fixed << setprecision(1) << t_now << "s)..." << endl;
        }
    }
    */
    
    // 输出轨迹
    for (size_t i = 0; i < full_traj.size(); ++i) {
        const auto& s = full_traj[i];
        Vector3d att = INSMath::m2att(INSMath::q2mat(s.qnb)) / glv.deg;
        Vector3d pos = s.pos; 
        pos(0) /= glv.deg; 
        pos(1) /= glv.deg;
        
        double t = (nav_start_idx + i) * ts;
        f_traj << setprecision(10) << t << "," 
               << pos(0) << "," << pos(1) << "," << pos(2) << ","
               << att(0) << "," << att(1) << "," << att(2) << endl;
    }
    
    f_bias.close();
    f_traj.close();
    
    // 统计结果（直接用 engine.ins 而不是 full_traj）
    {
        Vector3d p1;
        p1(0) = engine.ins.pos(0) / glv.deg;
        p1(1) = engine.ins.pos(1) / glv.deg;
        p1(2) = engine.ins.pos(2);
        
        double drift = (p1.head<2>() - p0.head<2>()).norm() * 111320.0;
        double actual_time_h = num_segments * SEGMENT_SIZE * ts / 3600.0;
        
        cout << "------------------------------------------------" << endl;
        cout << fixed << setprecision(6);
        cout << "Init Pos (deg): " << p0.transpose() << endl;
        cout << "End Pos (deg):  " << p1.transpose() << endl;
        cout << "Total Horizontal Drift: " << drift << " m" << endl;
        if (actual_time_h > 0)
            cout << "Drift Rate: " << (drift / 1852.0) / actual_time_h << " nm/h" << endl;
        cout << "Final Gyro Bias (deg/h): " << (current_state.eb / glv.deg * 3600.0).transpose() << endl;
    }

    return 0;
}