# Hybrid INS System - 多速率惯性导航系统仿真平台

## 📋 项目概述

本项目实现了一个**基于原子传感器辅助的多速率惯性导航系统（Hybrid INS）**仿真平台，用于研究和验证多速率传感器融合算法在长航时导航中的性能。

### 核心特性

- **多速率传感器融合**：融合高频光纤陀螺（FOG，400Hz）与低频原子陀螺（CAI，0.5Hz）
- **三种算法实现对比**：
  - 纯光纤惯导基线（Pure FOG INS）
  - 因子图解析平滑（ES-FGO，Event-Sparse Factor Graph Optimization）
  - 传统延迟状态EKF（Delayed-State EKF）
- **长航时仿真能力**：支持24小时连续导航仿真
- **完整的对准流程**：粗对准 + KF精对准 + 几何均值滤波
- **原子传感器仿真器**：高保真原子陀螺和原子加计物理模型

---

## 🏗️ 架构设计

### 目录结构

```
hybrid_ins_cpp/
├── src/                    # 源代码目录
│   ├── support.h          # 工具函数、数据结构、对准算法
│   ├── cai_sim.h/cpp      # 原子传感器仿真器
│   ├── test.cpp           # 批量实验脚本（死区时间研究）
│   ├── test_fgo.cpp       # 论文核心对比实验（三种算法）
│   └── main.cpp           # 早期版本主程序
│
├── psins/                  # 捷联惯导库（PSINS）
│   ├── glv.h              # 全局常量定义
│   ├── earth.h/cpp        # 地球模型（重力、地球自转等）
│   ├── ins_math.h/cpp     # 惯导数学函数（旋转变换等）
│   ├── ins_state.h/cpp    # 惯导状态（姿态、速度、位置）
│   ├── kf_state.h/cpp     # 卡尔曼滤波状态
│   └── sins_engine.h/cpp  # SINS引擎（导航解算、对准）
│
├── doc/                    # 文档目录
│   ├── README.md          # 本文档
│   ├── atom_gyro_smoothing.md
│   ├── fgo_experiments.md
│   └── optimal_params.md
│
├── build/                  # 编译输出目录
├── figs_paper_final/      # 论文图表
└── *.csv                   # 数据文件（被.gitignore过滤）
```

---

## 🔬 核心模块详解

### 1. 原子传感器仿真器（`cai_sim.h/cpp`）

#### 原子陀螺仿真器 (`AtomicGyroSimulator`)

**物理模型**：
- 测量地球自转在体坐标系的投影：$\omega_{ie}^b = C_n^b \omega_{ie}^n$
- 包含散粒噪声（角度随机游走 ARW）
- 包含零偏不稳定性

**关键参数**：
```cpp
struct CAIParams {
    double T_cycle = 2.0;        // 干涉周期 2秒
    double arw_dpsh = 0.0e-5;    // ARW (deg/√h)
    double bias_dph = 1.0e-5;    // 零偏不稳定性 (deg/h)
    int seed = 123;              // 随机种子
};
```

**使用方式**：
```cpp
// 1. 构造时传入位置
AtomicGyroSimulator atom(pos_rad, glv);

// 2. 用对准后的姿态初始化（计算固定的 wie_b）
atom.Init(att_align);

// 3. 获取测量值
Vector3d wb_atom = atom.Measure();  // rad/s
```

#### 原子加计仿真器 (`AtomicAccSimulator`)

**物理模型**：
- 测量比力（Specific Force）
- 使用锁定姿态计算理想参考比力：$f^b_{ref} = (C_n^b)^T \cdot (-g^n)$
- 误差模型：零偏 + 白噪声

**关键特性**：
- 初始化时锁定对准姿态，提供稳定的重力参考
- 可用于抑制舒拉振荡（Schuler Oscillation）

---

### 2. 混合对准算法（`support.h` + `sins_engine.cpp`）

#### 对准流程

```
[输入] 10分钟静态IMU数据
   ↓
[粗对准] 解析式双矢量定姿（60秒）
   ├─ 陀螺积分获取地球自转方向
   └─ 加计平均获取重力方向
   ↓
[KF精对准] 卡尔曼滤波（剩余时间）
   ├─ 状态：姿态误差角、速度、陀螺零偏、加计零偏
   ├─ 观测：速度约束（静基座，速度应为零）
   └─ 输出：精对准姿态和零偏初值
   ↓
[几何均值滤波] 零偏稳定性约束
   ├─ 陀螺：eb_final = eb_raw * scale_factor
   ├─ 加计：db_final = db_raw * scale_factor
   └─ scale = sqrt(sigma_allan / sigma_measured)
   ↓
[输出] 高精度初始姿态和零偏
```

**配置参数**：
```cpp
struct HybridAlignConfig {
    double t_coarse = 60.0;        // 粗对准时长（秒）
    double t_fine = 3600.0;        // 精对准时长（秒）
    double eb_sigma_allan = 0.003; // 陀螺Allan零偏稳定性 (deg/h)
    double db_sigma_allan = 50.0;  // 加计Allan零偏稳定性 (ug)
    bool verbose = true;
};
```

---

### 3. 三种导航算法实现（`test_fgo.cpp`）

#### 算法1：纯光纤惯导基线（Pure FOG INS）

**目的**：建立性能下界，测量无原子传感器约束时的最大漂移

**实现**：
```cpp
void RunPureFOG(const ExperimentConfig& exp, const SimulationContext& ctx) {
    // 初始化SINS引擎
    SinsEngine sinsegine(ctx.ts);
    sinsegine.ins = INSState(att_align, vel_zero, pos_ref, ...);
    
    // 循环：纯机械编排
    while (nav_loader.Next(epoch)) {
        sinsegine.Step_Nav(epoch.wm, epoch.vm);  // 惯性解算
        sinsegine.ins.vn(2) = 0.0;               // 高度阻尼
        
        // 记录漂移
        if (每1秒) {
            计算水平漂移距离
            保存到CSV
        }
    }
}
```

**特点**：
- 无外部观测约束
- 24小时漂移可达 **16.4 km**（典型值）
- 计算复杂度：O(1)
- 输出文件：`nav_pure_fog.csv`

---

#### 算法2：因子图解析平滑（ES-FGO，本文核心贡献）

**目的**：多速率融合，实时输出锯齿状轨迹，测算O(1)平滑算力

**核心思想**：
1. **周期性零偏估计**（每2秒）：
   ```
   有效时间 t_active = 1.6s → 累积FOG数据
   死区时间 t_dead = 0.4s   → 原子陀螺测量
   
   零偏估计 = FOG_均值 - 原子_测量
   ```

2. **雅可比累积**（实时）：
   ```cpp
   // 每个采样周期
   Phi_xb = [-Cnb*dt, 0, 0; 0, -Cnb*dt, 0; 0, 0, 0]
   Phi_xx = [I-[winn]*dt, 0, 0; [fn]*dt, I-[wcor]*dt, 0; 0, Mpv*dt, I]
   J = Phi_xx * J + Phi_xb  // 累积雅可比
   ```

3. **解析平滑更新**（周期末）：
   ```cpp
   // 计算零偏增量
   delta_b = [eb_new - eb_old; db_new - db_old]
   
   // O(1) 解析回溯历史状态
   for (历史时刻 t in log_buffer) {
       delta_x = J_at_t * delta_b
       姿态修正、速度修正、位置修正
   }
   
   // 更新当前时刻
   delta_x_now = J * delta_b
   更新姿态、速度、位置、零偏
   ```

**关键优势**：
- **O(1) 复杂度**：无论回溯多少历史时刻，计算量恒定
- **实时锯齿轨迹**：高频记录（10Hz）前向漂移和垂直拉回
- **零方差更新**：解析式平滑，无需维护协方差矩阵

**性能指标**：
- 24小时漂移：**< 50 m**（典型值）
- 平均更新耗时：**0.003 ms/cycle**
- 输出文件：`nav_esfgo_full.csv`

---

#### 算法3：传统延迟状态EKF（Delayed-State EKF）

**目的**：对比算力开销，展示O(N)重积分的计算负担

**实现**：
```cpp
void RunEKF(const ExperimentConfig& exp, const SimulationContext& ctx) {
    // 保存周期起始状态和原始IMU
    INSState state_t0 = sinsegine.ins;
    vector<IMUData> raw_imu_buffer;
    
    while (nav_loader.Next(epoch)) {
        raw_imu_buffer.push_back(epoch);
        标称前向传播...
        
        if (周期结束) {
            估计新零偏
            
            // --- 重积分开始 ---
            sinsegine.ins = state_t0;         // 回滚状态
            sinsegine.ins.set_bias(eb_new, db_new);
            
            for (ep in raw_imu_buffer) {      // O(N)
                sinsegine.Step_Nav(ep.wm, ep.vm);
            }
            // --- 重积分结束 ---
            
            更新 state_t0
            清空 raw_imu_buffer
        }
    }
}
```

**特点**：
- 计算复杂度：**O(N)**，N为周期内采样点数
- 1小时平均更新耗时：**0.126 ms/cycle**（约为ES-FGO的42倍）
- 无文件输出（仅用于算力测试）

---

### 4. 批量实验脚本（`test.cpp`）

用于系统性研究不同参数配置的影响，主要聚焦**死区时间（t_dead）**和**原子加计使能**对导航性能的影响。

**实验配置示例**：
```cpp
vector<ExperimentConfig> experiments = {
    // 仅原子陀螺（复现舒拉振荡）
    {"Group_1.6s_GyroOnly", ES_FGO, 24h, 1.6s, 0.4s, false, "...csv"},
    
    // 原子陀螺+加计（压制舒拉振荡）
    {"Group_1.6s_GyroAcc",  ES_FGO, 24h, 1.6s, 0.4s, true,  "...csv"},
};
```

**输出数据**（每秒记录）：
- 位置误差（纬度、经度、高度）
- 水平漂移距离
- 速度、姿态
- 零偏估计值
- 死区积分残差（用于验证算法正确性）

---

## 🔧 编译与运行

### 依赖项

- **C++11** 或更高版本
- **Eigen3** 线性代数库
- **CMake** 3.10+
- **Make** 或 Ninja

### 编译步骤

```bash
cd /home/v/dev/hybrid_ins_cpp
mkdir -p build && cd build
cmake ..
make -j4
```

### 运行程序

#### 1. 论文核心对比实验（推荐）

```bash
cd build
./cai_sim_test
```

**功能**：
- 自动执行混合对准（10分钟）
- 运行3个对比实验：
  1. Pure FOG（24小时）
  2. ES-FGO（24小时）
  3. EKF（1小时，仅测算力）
- 输出CSV文件到当前目录

**输出文件**：
- `nav_pure_fog.csv` - 纯光纤基线轨迹
- `nav_esfgo_full.csv` - ES-FGO完整轨迹（含锯齿）
- `nav_ekf_perf.csv` - EKF性能测试
- `nav_dead_1.6s_gyro_only.csv` - 消融实验（无加计）

---

#### 2. 批量参数研究

```bash
cd build
./test
```

**功能**：
- 对准 + 批量运行多组实验
- 可自定义死区时间、有效时间等参数
- 适合做参数扫描和灵敏度分析

**修改实验配置**：
编辑 `src/test.cpp` 第 221-224 行：
```cpp
vector<ExperimentConfig> experiments = {
    {"CustomExp1", ES_FGO, 24.0, 1.5, 0.5, true, "output1.csv"},
    {"CustomExp2", ES_FGO, 24.0, 1.8, 0.2, true, "output2.csv"},
    // ... 添加更多配置
};
```

---

## 📊 数据处理与可视化

### Python脚本（位于 `/home/v/dev/scripts/`）

#### 1. 绘制Allan方差图
```bash
cd /home/v/dev/scripts
# 陀螺和加计Allan方差（合并版）
python3 plot_allan_variance.py              # 默认绘制两者
python3 plot_allan_variance.py -t gyro      # 仅陀螺
python3 plot_allan_variance.py -t acc       # 仅加计
```

#### 2. 绘制导航结果
```bash
python3 plot_algorithm_comparison.py  # 三算法对比（Pure FOG vs ES-FGO）
python3 plot_schuler_oscillation.py   # 舒拉振荡分析（消融实验）
python3 plot_paper.py                 # 论文主图（轨迹、误差曲线）
python3 plot_dead_time.py             # 死区时间分析
```

### 输出图表位置
- `figs_paper_final/` - 论文最终图表
- `nav_dead_1.6s_gyro_acc_figs/` - 加计辅助实验图
- `nav_dead_1.6s_gyro_only_figs/` - 纯陀螺实验图

---

## 📐 关键物理模型

### 1. 地球模型（`earth.h/cpp`）

**重力模型**（GRS80椭球）：
```cpp
double sl = sin(lat);
double gL = ge * (1.0 + beta * sl * sl + beta1 * sl * sl * sl * sl);
gn = Vector3d(0, 0, -gL);  // 指向地心
```

**地球自转**：
```cpp
wien = Vector3d(0, wie * cos(lat), wie * sin(lat));  // n系
winn = Cnb * wien;  // b系投影
```

**哥氏加速度**：
```cpp
wcor = 2 * wien + wen;  // wen为牵连角速度
```

### 2. 姿态更新（四元数法）

```cpp
// 等效旋转矢量
Vector3d phi = wm - (eb + wien*dt) * dt;

// 四元数更新
Quaterniond dq = INSMath::rv2q(phi);
qnb = qnb * dq;
qnb.normalize();

// 转姿态矩阵
Cnb = INSMath::q2mat(qnb);
att = INSMath::m2att(Cnb);
```

### 3. 速度更新

```cpp
Vector3d fn = Cnb * (vm - db*dt);  // 比力转n系
Vector3d dvn = fn - (2*wien + wen).cross(vn) + gn*dt;
vn += dvn;
vn(2) = 0;  // 高度阻尼（可选）
```

### 4. 位置更新

```cpp
Matrix3d Mpv;
Mpv(0,1) = 1.0 / (RM + h);              // 纬度
Mpv(1,0) = 1.0 / ((RN + h) * cos(lat)); // 经度
Mpv(2,2) = 1.0;                         // 高度

pos += Mpv * vn * dt;
```

---

## 🎯 实验场景配置

### 参考位置
- **纬度**：32.0286° N
- **经度**：118.8533° E
- **高度**：17.0 m
- **地理位置**：江苏省南京市

### 传感器参数

#### 光纤陀螺（FOG）
- 采样率：400 Hz
- 零偏不稳定性：0.003 deg/h（Allan方差）
- 角度随机游走：量级 deg/√h

#### 光纤加计
- 采样率：400 Hz
- 零偏不稳定性：50 ug（Allan方差）
- 速度随机游走：量级 ug

#### 原子陀螺（CAI）
- 更新周期：2.0 s
- 零偏不稳定性：1e-5 deg/h
- 角度随机游走：接近零

#### 原子加计
- 更新周期：2.0 s
- 零偏：可配置（通常设为0）
- 白噪声：0.05 ug

### 多速率配置
- **有效时间** t_active = 1.6 s
- **死区时间** t_dead = 0.4 s
- **总周期** T_cycle = 2.0 s
- **FOG采样点数**：800点/周期

---

## 📈 性能指标对比

| 算法 | 24h漂移 | 更新耗时 | 计算复杂度 | 历史平滑 |
|------|---------|----------|-----------|----------|
| Pure FOG | 16.4 km | - | O(1) | 无 |
| ES-FGO | < 50 m | 0.003 ms | O(1) | 解析式 |
| Delayed-State EKF | < 50 m | 0.126 ms | O(N) | 重积分 |

**ES-FGO优势**：
- ✅ 零方差平滑（无需协方差传播）
- ✅ 算力极低（是EKF的1/42）
- ✅ 实时输出（适合在线应用）
- ✅ 可扩展性强（历史窗口可任意长）

---

## 🔍 故障排查

### 常见问题

**Q1: 编译报错 "cannot find -lEigen3"**

```bash
# Ubuntu/Debian
sudo apt-get install libeigen3-dev

# 或手动指定路径
cmake -DEIGEN3_INCLUDE_DIR=/usr/include/eigen3 ..
```

**Q2: 运行时提示 "Cannot open fog_part1.csv"**

CSV数据文件位于项目根目录，运行程序需在 `build/` 目录下执行：
```bash
cd build
./cai_sim_test  # 正确
```

**Q3: 对准失败 "Alignment failed!"**

- 检查数据文件完整性（应有完整6小时静态数据）
- 调整对准配置参数（`HybridAlignConfig`）
- 增加verbose输出查看中间结果

**Q4: 导航结果异常发散**

- 确认对准成功（输出显示姿态和零偏）
- 检查原子传感器初始化（必须用对准姿态Init）
- 验证数据时间戳连续性

---

---

## 🤝 贡献指南

### 代码规范
- 使用4空格缩进
- 变量命名：蛇形命名法（`local_g`）
- 类名：大驼峰（`AtomicGyroSimulator`）
- 注释：中英文混合，关键算法必须详细注释

### 提交信息格式
```
类型(作用域): 简短描述

详细说明（可选）

关联Issue: #123
```

类型：feat（新功能）、fix（修复）、docs（文档）、perf（性能）

---

## 📝 引用

如果本项目对您的研究有帮助，请引用：

```bibtex
@article{YourPaper2026,
  title={Multi-Rate Inertial Navigation System with Atomic Sensors},
  author={Your Name},
  journal={IEEE Transactions on ...},
  year={2026}
}
```

---

## 📄 许可证

本项目遵循 MIT 许可证。详见 [LICENSE](../LICENSE) 文件。

---

## 联系方式

- **作者**：v
- **邮箱**：v@example.com
- **项目主页**：https://github.com/flowercn/hybird_ins_cpp

---

**最后更新**：2026年3月4日
