# 因子图优化 (FGO) 实验总结

## 1. 概述

本文档记录了使用 GTSAM 因子图优化改进惯导精度的实验过程，包括成功案例和失败教训。

## 2. 实验环境

- **数据源**: fog3h.csv (3小时静态FOG数据，400Hz)
- **位置**: 南京 (32.0286°N, 118.8533°E)
- **FOG参数**: ARW = 0.003 deg/√h, 刻度因数 K = 1.000372 (371.67 ppm)
- **原子陀螺参数**: ARW = 2e-4 deg/√h, T_cycle = 2s, Bias = 1e-5 deg/h

## 3. 基线性能

### 3.1 传统纯惯导
- **漂移率**: 1.08 nm/h (使用最优对准参数)
- **方法**: KF精对准 + 几何均值零偏估计

### 3.2 最优对准参数 (Grid Search获得)
```
姿态 (deg): [0.0344, 0.335, 0.590]
零偏 (deg/h): [-0.0017, -0.0026, -0.0103]
```

## 4. 成功实验: FOG-only FGO

### 4.1 实现思路
```
因子图结构:
  X0 ---- preint ---- X1 ---- preint ---- X2 ...
  |                   |                   |
  prior               zero_vel            zero_vel
                      zero_pos            zero_pos
```

- **节点间隔**: 1秒
- **约束**: 零速约束 + 零位移约束 (静态场景)
- **预积分**: CombinedImuFactor

### 4.2 关键参数
```cpp
// IMU噪声参数 (连续时间)
gyro_noise_density = 0.003 * deg2rad / 60.0;   // rad/s/√Hz
accel_noise_density = 50e-6 * 9.8;              // m/s²/√Hz
gyro_bias_rw = 0.0003 * deg2rad / 3600.0;       // rad/s²/√Hz
accel_bias_rw = 1e-6 * 9.8;                     // m/s³/√Hz

// 约束噪声
zero_vel_sigma = 0.001;  // m/s
zero_pos_sigma = 0.01;   // m
```

### 4.3 结果
| 方法 | 漂移率 | 相对提升 |
|------|--------|----------|
| 传统INS | 1.08 nm/h | 基线 |
| **FGO (FOG-only)** | **0.72 nm/h** | **33.6%** |

### 4.4 成功原因
1. 零速/零位移约束提供了强观测信息
2. 批量优化可以利用全局信息
3. GTSAM 的预积分正确处理了IMU误差传播

## 5. 失败实验: 混合FGO (FOG + 原子陀螺)

### 5.1 尝试1: fgo_hybrid_test.cpp
**思路**: 添加地球自转约束，用原子陀螺测量的 wie_b 作为参考

```cpp
// 约束: FOG测量的角速度应等于地球自转
// wie_b_ref = Cnb^T * wie_n (从对准姿态计算)
BetweenFactor<Vector3>(key, wie_b_ref, noise)
```

**结果**: 491 nm/h (灾难性失败)

**失败原因**:
- wie_b_ref 是基于"理想"对准姿态计算的
- FOG实际测量 = wie_b + eb (含零偏)
- 两者不一致导致优化器试图"修正"姿态来匹配，反而引入巨大误差

### 5.2 尝试2: fgo_hybrid_v2.cpp
**思路**: 调整权重，弱化地球自转约束

```cpp
// 降低约束权重
earth_rotation_sigma = 0.01 * deg2rad / 3600.0;  // 很大的不确定性
```

**结果**: 4.5 nm/h (比传统方法还差)

**失败原因**:
- 约束太弱等于没用
- 约束太强又会干扰优化
- 本质问题没解决

### 5.3 尝试3: fgo_hybrid_v3.cpp
**思路**: 只用原子陀螺约束估计零偏，不约束姿态

```cpp
// 零偏先验因子
// eb_prior = FOG_avg - wie_b_ref
PriorFactor<imuBias>(B(0), eb_prior, noise)
```

**结果**: 144 nm/h (严重恶化)

**失败原因**:
- eb_prior 的计算依然依赖 wie_b_ref
- 而 wie_b_ref 的精度取决于对准姿态
- 形成循环依赖

### 5.4 尝试4: fgo_atom_test.cpp
**思路**: 激进方案，每个节点都添加原子陀螺测量约束

```cpp
// 每秒都约束角速度
for (each node) {
    graph.add(PriorFactor(omega, wie_b_ref + noise))
}
```

**结果**: 8.7×10⁹ nm/h (完全发散)

**失败原因**:
- 过约束导致优化器无法收敛
- 约束之间相互矛盾

## 6. 核心教训

### 6.1 根本问题
**原子陀螺约束的核心矛盾**:
1. 原子陀螺测量的是"真实"地球自转 wie_b = Cnb_true^T * wie_n
2. FOG测量的是 wm = wie_b + eb (含零偏和噪声)
3. 要建立约束，需要知道 Cnb_true
4. 但 Cnb_true 正是我们想要估计的！

### 6.2 失败模式
```
                    ┌─────────────┐
         ┌─────────►│ 对准姿态估计 │
         │          └──────┬──────┘
         │                 │
         │                 ▼ (用于计算)
         │          ┌─────────────┐
         │          │   wie_b_ref  │
         │          └──────┬──────┘
         │                 │
         │                 ▼ (作为约束)
         │          ┌─────────────┐
         └──────────│  FGO 优化   │◄─── 循环依赖!
                    └─────────────┘
```

### 6.3 为什么回溯平滑可行而FGO不行？

| 方面 | 回溯平滑 | FGO约束 |
|------|----------|---------|
| wie_b来源 | 已知真值 (仿真) | 从对准估计 |
| 误差处理 | 差分消除 | 绝对值约束 |
| 容错性 | 高 | 低 |

回溯平滑用的是 `eb = FOG_avg - wie_b_true`，如果 wie_b_true 不准，eb 估计也不准，但至少方向一致。

FGO约束用的是 `FOG应该等于wie_b_ref`，如果 wie_b_ref 不准，优化器会扭曲整个状态来满足约束。

## 7. 未来方向

### 7.1 可能可行的方案
1. **松耦合**: 原子陀螺只用于后处理校正，不进入FGO
2. **增量约束**: 不约束绝对值，只约束变化量 (类似IMU预积分)
3. **联合估计**: 把 wie_b 也作为状态量同时估计

### 7.2 需要的数学推导
- 原子陀螺观测方程的正确形式
- 与IMU状态的耦合方式
- 可观测性分析

## 8. 项目文件结构

### 8.1 核心代码 (src/)
| 文件 | 说明 |
|------|------|
| `cai_sim.h` | 原子陀螺仿真器头文件 |
| `cai_sim.cpp` | 原子陀螺仿真器实现 |
| `eskf_cai.h` | ESKF 15维滤波器 |
| `support.h` | 数据加载、辅助函数 |
| `test.cpp` | 主测试程序 |
| `hybirdmain.cpp` | 混合导航主程序 (开发中) |

### 8.2 PSINS库 (psins/)
| 文件 | 说明 |
|------|------|
| `glv.h` | 全局常量 (地球参数、单位转换) |
| `earth.h/cpp` | 地球模型 (重力、曲率半径) |
| `ins_math.h/cpp` | 数学工具 (四元数、旋转矩阵) |
| `ins_state.h/cpp` | INS状态定义 |
| `kf_state.h/cpp` | 卡尔曼滤波器 |
| `sins_engine.h/cpp` | 惯导解算引擎 |

### 8.3 归档代码 (archive/)
| 文件 | 说明 |
|------|------|
| `fgo_test.cpp` | ⭐ FGO成功案例 (0.72 nm/h) |
| `fgo_sim_test.cpp` | FGO仿真验证 |
| `smooth_demo.cpp` | 回溯平滑演示 |
| `strategy_compare.cpp` | 策略对比 |
| `theory_test.cpp` | 理论极限分析 |
| `example_minimal.cpp` | 最简使用示例 |
| `example_complete.cpp` | 完整使用示例 |

### 8.4 已删除的失败实验 (不保留)
| 文件 | 结果 | 删除原因 |
|------|------|----------|
| `fgo_hybrid_test.cpp` | 491 nm/h | 地球自转约束错误 |
| `fgo_hybrid_v2.cpp` | 4.5 nm/h | 权重调整无效 |
| `fgo_hybrid_v3.cpp` | 144 nm/h | 零偏先验错误 |
| `fgo_atom_test.cpp` | 发散 | 过约束 |
| `fgo_core.h/cpp` | - | 未完成的封装 |
| `fgo_main.cpp` | - | 废弃入口 |

### 8.5 数据文件
| 文件 | 说明 |
|------|------|
| `fog3h.csv` | 3小时FOG静态数据 (400Hz) |
| `fog.csv` | 1小时FOG数据 |

### 8.6 文档 (doc/)
| 文件 | 说明 |
|------|------|
| `optimal_params.md` | 最优对准参数记录 |
| `fgo_experiments.md` | 本文档 |
| `atom_gyro_smoothing.md` | 原子陀螺平滑总结 |

## 9. 代码示例

成功的 FGO 核心代码片段 (来自 `archive/fgo_test.cpp`):

```cpp
// 预积分参数
auto preint_params = PreintegrationCombinedParams::MakeSharedU(local_g);
preint_params->setGyroscopeCovariance(pow(gyro_noise_density, 2) * Matrix3d::Identity());
preint_params->setAccelerometerCovariance(pow(accel_noise_density, 2) * Matrix3d::Identity());

// 添加因子
graph.emplace_shared<CombinedImuFactor>(
    X(i), V(i), X(i+1), V(i+1), B(i), B(i+1), *preint);

// 零速约束
graph.addPrior(V(i+1), zero_vel, noise_zero_vel);

// 零位移约束 (只约束位置，不约束姿态)
auto noise_pos_only = noiseModel::Diagonal::Sigmas(
    (Vector6() << 1e3, 1e3, 1e3, 0.01, 0.01, 0.01).finished());
graph.addPrior(X(i+1), pose0, noise_pos_only);
```

## 10. 参考文献

- GTSAM Tutorial: IMU Preintegration
- Forster et al., "On-Manifold Preintegration for Real-Time Visual-Inertial Odometry"
