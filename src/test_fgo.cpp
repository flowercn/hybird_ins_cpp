#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include "sins_engine.h"
#include "cai_sim.h"
#include "support.h" 

using namespace std;
using namespace Eigen;
using namespace cai;

// 辅助函数：生成反对称矩阵
auto askew = [](const Vector3d& v) -> Matrix3d {
    Matrix3d m;
    m <<   0,  -v(2),  v(1),
         v(2),    0,  -v(0),
        -v(1),  v(0),    0;
    return m;
};

// ==========================================
// 用于 FGO 解析回溯的日志缓存结构
// ==========================================
struct LogPoint {
    double nav_time;
    // 直接缓存基础类型，彻底避开 INSState 的构造函数限制
    Vector3d att;
    Vector3d vn;
    Vector3d pos;
    Matrix3d Cnb;
    MatrixXd J_at_t;
    
    // 提供基础构造函数
    LogPoint(double t, const Vector3d& a, const Vector3d& v, const Vector3d& p, const Matrix3d& C, const MatrixXd& j) 
        : nav_time(t), att(a), vn(v), pos(p), Cnb(C), J_at_t(j) {}
};

// ==========================================
// 单步运行函数
// ==========================================
void RunSingleExperiment(const ExperimentConfig& exp, const SimulationContext& ctx) {
    cout << "\n--------------------------------------------------" << endl;
    cout << " Running Scenario: " << exp.name << endl;
    cout << " Config: Active=" << exp.t_active << "s, Dead=" << exp.t_dead << "s -> " << exp.output_file << endl;
    cout << "--------------------------------------------------" << endl;

    // A. 初始化传感器
    AtomicGyroSimulator atom(ctx.pos_ref, ctx.glv);
    atom.Init(ctx.att_align); 
    
    CAIParams acc_params;
    acc_params.bias_ug = 0.0; 
    acc_params.vrw_ug  = 0.05; 
    AtomicAccSimulator atom_acc(ctx.pos_ref, ctx.glv, acc_params);
    atom_acc.Init(ctx.att_align); 
    
    // B. 初始化惯导引擎
    SinsEngine sinsegine(ctx.ts);
    sinsegine.eth = ctx.eth; 
    sinsegine.ins = INSState(ctx.att_align, Vector3d::Zero(), ctx.pos_ref, ctx.ts, sinsegine.eth);
    sinsegine.ins.set_bias(ctx.eb_align, ctx.db_align);

    std::vector<std::string> rest_files(ctx.file_list.begin() + 1, ctx.file_list.end());
    IMUChainedLoader nav_loader(rest_files, ctx.ts, ctx.local_g);

    ofstream log(exp.output_file);
    log << "time,lat_err,lon_err,h_err,drift,vn,ve,vu,roll,pitch,yaw,"
        << "eb_x,eb_y,eb_z,db_x,db_y,db_z,"
        << "res_gyro_x,res_gyro_y,res_gyro_z,"
        << "res_acc_x,res_acc_y,res_acc_z\n"; 

    double T_cycle = 2.0; 
    int total_samples_per_cycle = static_cast<int>(T_cycle / ctx.ts); 
    int active_samples = static_cast<int>(exp.t_active / ctx.ts);     
    
    Vector3d fog_gyro_active_sum = Vector3d::Zero(); 
    Vector3d fog_acc_active_sum  = Vector3d::Zero(); 
    Vector3d gyro_dead_raw_sum = Vector3d::Zero(); 
    Vector3d acc_dead_raw_sum  = Vector3d::Zero(); 
    
    int sample_count = 0;
    int dead_count = 0;
    double max_drift = 0;
    size_t total_steps = 0;

    // ==========================================
    // FGO 核心变量：超级雅可比矩阵与日志缓存
    // ==========================================
    MatrixXd J = MatrixXd::Zero(9, 6); 
    MatrixXd Phi_xx = MatrixXd::Identity(9, 9);
    MatrixXd Phi_xb = MatrixXd::Zero(9, 6);
    std::vector<LogPoint> log_buffer;
    log_buffer.reserve(10); 

    auto ProcessEpoch = [&](const IMUData& epoch) {
        sample_count++;
        total_steps++;

        // 1. 数据累加
        if (sample_count <= active_samples) {
            fog_gyro_active_sum += epoch.wm / ctx.ts; 
            fog_acc_active_sum  += epoch.vm / ctx.ts; 
        } else {
            Vector3d earth_inc_body = ctx.Cnb_align.transpose() * ctx.eth.wien * ctx.ts;
            gyro_dead_raw_sum += (epoch.wm - earth_inc_body);
            Vector3d force_inc_body = ctx.Cnb_align.transpose() * (-ctx.eth.gn) * ctx.ts;
            acc_dead_raw_sum += (epoch.vm - force_inc_body);
            dead_count++;
        }

        // 2. 标称物理状态向前推演 (绝不回头)
        sinsegine.Step_Nav(epoch.wm, epoch.vm);
        sinsegine.ins.vn(2) = 0.0; // 阻尼

        // 3. 伴随计算 FGO 雅可比矩阵
        Matrix3d Cnb = sinsegine.ins.Cnb;
        Vector3d fn = Cnb * (epoch.vm / ctx.ts);
        Vector3d winn = sinsegine.eth.winn;
        Vector3d wcor = sinsegine.eth.wcor;

        double RM_h = sinsegine.eth.RMh; // 统一使用地球参数的曲率半径
        double RN_h = sinsegine.eth.RNh;
        double cos_lat = cos(sinsegine.ins.pos(0));
        Matrix3d Mpv = Matrix3d::Zero();
        Mpv(0, 1) = 1.0 / RM_h;
        Mpv(1, 0) = 1.0 / (RN_h * cos_lat);
        Mpv(2, 2) = 1.0;

        // 驱动矩阵 (误差注入)
        Phi_xb.setZero();
        Phi_xb.block<3,3>(0, 0) = -Cnb * ctx.ts;
        Phi_xb.block<3,3>(3, 3) = -Cnb * ctx.ts;

        // 演化矩阵 (舒拉基因与地球自转)
        Phi_xx.setIdentity();
        Phi_xx.block<3,3>(0, 0) -= askew(winn) * ctx.ts;
        Phi_xx.block<3,3>(3, 0) += askew(fn) * ctx.ts;
        Phi_xx.block<3,3>(3, 3) -= askew(wcor) * ctx.ts;
        Phi_xx.block<3,3>(6, 3) += Mpv * ctx.ts;

        // 连乘累加雅可比
        J = Phi_xx * J + Phi_xb;

        // 4. 记录输出锚点 (提取基础数据类型存入 Buffer)
        if (total_steps % 400 == 0) { 
            log_buffer.emplace_back(
                total_steps * ctx.ts, 
                sinsegine.ins.att, 
                sinsegine.ins.vn, 
                sinsegine.ins.pos, 
                sinsegine.ins.Cnb, 
                J
            );
        }

        // 5. 周期结算与 FGO 解析平滑
        if (sample_count >= total_samples_per_cycle) {
            
            // 估算原子观测
            Vector3d fog_gyro_mean = Vector3d::Zero();
            Vector3d fog_acc_mean  = Vector3d::Zero();
            if (active_samples > 0) {
                fog_gyro_mean = fog_gyro_active_sum / active_samples;
                fog_acc_mean  = fog_acc_active_sum  / active_samples;
            }

            Vector3d atom_gyro_meas = atom.Measure(); 
            Vector3d eb_new = fog_gyro_mean - atom_gyro_meas;

            Vector3d db_new = sinsegine.ins.db;
            if (exp.use_atomic_acc && active_samples > 0) {
                Vector3d atom_acc_meas  = atom_acc.Measure(); 
                db_new = fog_acc_mean  - atom_acc_meas;
            }

            // 计算零偏真值与当前标称值的偏差向量 (Delta b)
            VectorXd delta_b(6);
            delta_b.segment<3>(0) = eb_new - sinsegine.ins.eb;
            delta_b.segment<3>(3) = db_new - sinsegine.ins.db;

            // 计算死区残差 (仅用于日志)
            double t_dead_real = dead_count * ctx.ts;
            Vector3d final_gyro_res = gyro_dead_raw_sum - eb_new * t_dead_real;
            Vector3d final_acc_res  = acc_dead_raw_sum  - db_new * t_dead_real;

            // --- 利用雅可比平滑历史输出，替代 for 回溯 ---
            for (const auto& lp : log_buffer) {
                VectorXd dx = lp.J_at_t * delta_b;
                
                // 严密姿态更新：利用 Eigen 的 AngleAxisd 执行精准李群映射
                Vector3d phi = dx.segment<3>(0);
                Matrix3d R_phi;
                if (phi.norm() > 1e-12) {
                    R_phi = AngleAxisd(phi.norm(), phi.normalized()).toRotationMatrix();
                } else {
                    R_phi = Matrix3d::Identity() + askew(phi); // 极小角退化保护
                }
                Matrix3d Cnb_cor = R_phi * lp.Cnb;     
                Vector3d s_att = INSMath::m2att(Cnb_cor);          
                
                Vector3d s_vn  = lp.vn  + dx.segment<3>(3);
                Vector3d s_pos = lp.pos + dx.segment<3>(6);

                double lat_err = (s_pos(0) - ctx.pos_ref(0)) * ctx.glv.Re;
                double lon_err = (s_pos(1) - ctx.pos_ref(1)) * ctx.glv.Re * cos(ctx.pos_ref(0));
                double h_err = s_pos(2) - ctx.pos_ref(2);
                double drift = sqrt(lat_err*lat_err + lon_err*lon_err);
                if(drift > max_drift) max_drift = drift;

                log << fixed << setprecision(15) 
                    << lp.nav_time << "," << lat_err << "," << lon_err << "," << h_err << "," << drift << ","
                    << s_vn(0) << "," << s_vn(1) << "," << s_vn(2) << "," 
                    << s_att(0)/ctx.glv.deg << "," << s_att(1)/ctx.glv.deg << "," << s_att(2)/ctx.glv.deg << ","
                    << eb_new(0)/ctx.glv.deg*3600 << "," << eb_new(1)/ctx.glv.deg*3600 << "," << eb_new(2)/ctx.glv.deg*3600 << ","
                    << db_new(0)/ctx.glv.ug << "," << db_new(1)/ctx.glv.ug << "," << db_new(2)/ctx.glv.ug << ","
                    << final_gyro_res(0) << "," << final_gyro_res(1) << "," << final_gyro_res(2) << ","
                    << final_acc_res(0)  << "," << final_acc_res(1)  << "," << final_acc_res(2) << "\n";
                
                if (lp.nav_time > 0 && static_cast<int>(lp.nav_time) % 3600 == 0) {
                     cout << "   Time: " << lp.nav_time / 3600.0 << " h | Drift: " << drift 
                          << " m | Max: " << max_drift << " m" << endl;
                }
            }

            // --- 修正 SinsEngine 当前状态，为下一个 2.0s 做准备 ---
            VectorXd final_dx = J * delta_b;
            
            // 当前姿态同步进行严密矩阵级修正
            Vector3d final_phi = final_dx.segment<3>(0);
            Matrix3d R_final_phi;
            if (final_phi.norm() > 1e-12) {
                R_final_phi = AngleAxisd(final_phi.norm(), final_phi.normalized()).toRotationMatrix();
            } else {
                R_final_phi = Matrix3d::Identity() + askew(final_phi);
            }
            Matrix3d Cnb_cor_final = R_final_phi * sinsegine.ins.Cnb;
            
            sinsegine.ins.att = INSMath::m2att(Cnb_cor_final);
            sinsegine.ins.Cnb = Cnb_cor_final;
            sinsegine.ins.qnb = INSMath::a2qua(sinsegine.ins.att); 
            
            sinsegine.ins.vn  += final_dx.segment<3>(3);
            sinsegine.ins.pos += final_dx.segment<3>(6);
            sinsegine.ins.set_bias(eb_new, db_new);

            // 清理与重置
            fog_gyro_active_sum.setZero(); fog_acc_active_sum.setZero();
            gyro_dead_raw_sum.setZero();   acc_dead_raw_sum.setZero();
            dead_count = 0;
            sample_count = 0; 
            J.setZero(); 
            log_buffer.clear();
        }
    };

    IMUData epoch;
    while (nav_loader.Next(epoch)) ProcessEpoch(epoch);
    log.close();
    cout << " [Done] Final Drift: " << max_drift << " m" << endl;
}

// ==========================================
// 主函数
// ==========================================
int main() {
    SimulationContext ctx;
    ctx.ts = 1.0 / 400.0;
    ctx.glv = GLV(); 
    ctx.pos_ref = Vector3d(32.0286 * ctx.glv.deg, 118.8533 * ctx.glv.deg, 17.0);
    ctx.eth = Earth(ctx.glv);
    ctx.eth.update(ctx.pos_ref, Vector3d::Zero());
    ctx.local_g = ctx.eth.gn.norm();
    
    ctx.file_list = {
        "../fog_part1.csv", "../fog_part2.csv", "../fog_part3.csv", 
        "../fog_part4.csv", "../fog_part5.csv"
    };

    cout << "\n[Phase 1] Loading Alignment Data: " << ctx.file_list[0] << endl;
    auto buffer_part1 = LoadIMUData(ctx.file_list[0], ctx.ts, ctx.local_g);
    
    double start_time_sec = 5.0 * 3600.0 + 50.0 * 60.0; 
    size_t start_idx = static_cast<size_t>(start_time_sec / ctx.ts);
    std::vector<IMUData> buffer_align(buffer_part1.begin() + start_idx, buffer_part1.end());

    SinsEngine align_engine(ctx.ts);
    align_engine.eth = ctx.eth;
    align_engine.ins = INSState(Vector3d::Zero(), Vector3d::Zero(), ctx.pos_ref, ctx.ts, ctx.eth);
    align_engine.res_init.pos = ctx.pos_ref;

    HybridAlignConfig cfg;
    cfg.t_coarse = 60;
    cfg.t_fine = buffer_align.size() * ctx.ts; 
    cfg.eb_sigma_allan = 0.003;
    cfg.db_sigma_allan = 50.0;
    cfg.verbose = false; 

    cout << "Running Alignment (" << cfg.t_fine << " s)..." << endl;
    auto align_res = align_engine.Run_HybridAlign(buffer_align, cfg);
    if (!align_res.valid) { cerr << "Alignment failed!" << endl; return -1; }
    
    cout << "=== ALIGNMENT SUCCESS ===" << endl;
    cout << "Att (deg):  " << align_res.att.transpose() / ctx.glv.deg << endl;
    cout << "Acc Bias :  " << align_res.db.transpose() / ctx.glv.ug << " ug" << endl;

    ctx.att_align = align_res.att;
    ctx.eb_align  = align_res.eb;
    ctx.db_align  = align_res.db;
    ctx.Cnb_align = INSMath::a2mat(align_res.att); 

    vector<IMUData>().swap(buffer_part1); 
    vector<IMUData>().swap(buffer_align);

    vector<ExperimentConfig> experiments = {
         {"Group_1.6s_GyroOnly", 1.6, 0.4, false, "nav_dead_1.6s_gyro_only.csv"}, 
         {"Group_1.6s_GyroAcc",  1.6, 0.4, true,  "nav_dead_1.6s_gyro_acc.csv"}, 
    };

    cout << "\n>>> Starting FGO Batch Processing (" << experiments.size() << " scenarios)..." << endl;
    
    for (const auto& exp : experiments) {
        RunSingleExperiment(exp, ctx);
    }

    cout << "\nAll Scenarios Completed Successfully!" << endl;
    return 0;
}