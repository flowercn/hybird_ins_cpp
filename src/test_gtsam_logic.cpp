#include <iostream>
#include <gtsam/geometry/Rot3.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/navigation/ImuFactor.h>
#include <gtsam/navigation/CombinedImuFactor.h>
#include <gtsam/base/Vector.h>

using namespace std;
using namespace gtsam;

// 辅助打印函数
void print_ypr(const string& label, const Vector3& ypr) {
    cout << label << ": [" 
         << ypr(0) * 180.0/M_PI << ", "   // Yaw
         << ypr(1) * 180.0/M_PI << ", "   // Pitch
         << ypr(2) * 180.0/M_PI << "] (deg)" << endl;
}

void test_rotation_order() {
    cout << "\n========== TEST 1: Rotation Order ==========" << endl;
    
    // 构造一个已知的欧拉角：Yaw=90, Pitch=30, Roll=0
    // GTSAM Rot3::Ypr(y, p, r) 
    double y_in = 90.0 * M_PI / 180.0;
    double p_in = 30.0 * M_PI / 180.0;
    double r_in = 10.0 * M_PI / 180.0;
    
    Rot3 R = Rot3::Ypr(y_in, p_in, r_in);
    
    // 获取 ypr() 的输出
    Vector3 ypr_out = R.ypr();
    
    print_ypr("Input (Expected)", Vector3(y_in, p_in, r_in));
    print_ypr("Output (Rot3::ypr)", ypr_out);
    
    if (std::abs(ypr_out(0) - y_in) < 1e-5 && 
        std::abs(ypr_out(1) - p_in) < 1e-5 && 
        std::abs(ypr_out(2) - r_in) < 1e-5) {
        cout << "✅ Conclusion: Rot3::ypr() returns [Yaw, Pitch, Roll] strictly." << endl;
    } else {
        cout << "❌ Conclusion: Order is different!" << endl;
    }
}

void test_integration_frame() {
    cout << "\n========== TEST 2: Integration Frame (Local vs ECEF) ==========" << endl;

    // 配置 IMU 参数
    // 假设我们在北向运动，重力向下 (NED系)
    // 如果是 ECEF，重力向量会随位置变化；如果是局部系，重力是常数
    Vector3 g_vec(0, 0, 9.81); // 假设 Z 轴向下 (NED) 的重力补偿
    
    auto params = PreintegratedCombinedMeasurements::Params::MakeSharedD(9.81);
    // 强制指定重力向量方向，看看积分出来的位移方向
    params->n_gravity = Vector3(0, 0, 9.81); 
    
    imuBias::ConstantBias bias; // 零偏为0
    PreintegratedCombinedMeasurements pim(params, bias);

    // 模拟 1秒 的数据
    double dt = 0.1;
    for(int i=0; i<10; ++i) {
        // IMU 读数：
        // 加速度计测得 [1.0, 0, -9.81] 
        // 解释：Z轴测得 -9.81 是因为支持力向上，或者单纯重力分量
        // 如果我们在 NED 系静止，Acc读数通常是 [0, 0, -9.8] (视传感器定义而定)
        // 这里我们要模拟：向 X 轴加速 1m/s^2，同时抵消重力
        // Acc Meas = a_body + R_expected * g_local
        // 假设 flat earth, R=I. Acc_meas = [1.0, 0, 0] - [0, 0, 9.81] = [1.0, 0, -9.81]
        
        Vector3 acc_meas(1.0, 0.0, -9.81); 
        Vector3 gyro_meas(0, 0, 0); // 无转动
        
        pim.integrateMeasurement(acc_meas, gyro_meas, dt);
    }

    // 预测结果
    // 初始状态：原点，静止
    Pose3 start_pose = Pose3();
    Vector3 start_vel = Vector3::Zero();
    
    NavState result = pim.predict(NavState(start_pose, start_vel), bias);
    
    cout << "Initial Pos: " << start_pose.translation().transpose() << endl;
    cout << "Final Pos:   " << result.pose().translation().transpose() << endl;
    
    // 理论计算 (局部笛卡尔系):
    // a_n = 1.0 m/s^2
    // t = 1.0 s
    // dist = 0.5 * a * t^2 = 0.5 * 1.0 * 1.0 = 0.5 meter
    
    if (std::abs(result.pose().x() - 0.5) < 0.05) {
        cout << "✅ Conclusion: GTSAM integrates in LOCAL CARTESIAN frame (Relative Displacement)." << endl;
        cout << "   Output is (x, y, z) in meters relative to Start Pose." << endl;
    } else {
        cout << "❌ Conclusion: Unexpected behavior. Result X=" << result.pose().x() << endl;
    }
}

int main() {
    test_rotation_order();
    test_integration_frame();
    return 0;
}