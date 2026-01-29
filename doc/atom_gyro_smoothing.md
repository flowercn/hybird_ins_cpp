# 原子陀螺回溯平滑实验总结

## 1. 概述

本文档记录了使用原子陀螺辅助FOG惯导的回溯平滑算法实验，包括算法原理、实验结果和参数分析。

## 2. 基本原理

### 2.1 核心思想
原子陀螺 (Cold Atom Interferometer, CAI) 能精确测量地球自转角速度，可作为FOG零偏估计的参考基准。

```
FOG测量:  wm_fog = wie_b + eb + noise_fog
原子测量: wm_cai = wie_b + noise_cai (无零偏)

差分:     eb ≈ wm_fog - wm_cai
```

### 2.2 回溯平滑流程
```
          第一遍 (收集)              第二遍 (回放)
        ─────────────────        ─────────────────
时间 →  |  T_cycle  |  T_cycle  |    重新导航
        |           |           |
        ├───────────┼───────────┤
        │ FOG积分   │ FOG积分   │    应用校正后的 eb
        │ 原子测量  │ 原子测量  │
        │ → eb_1    │ → eb_2    │
        └───────────┴───────────┘
```

## 3. 设备参数

### 3.1 FOG (光纤陀螺)
| 参数 | 值 | 说明 |
|------|-----|------|
| ARW | 0.003 deg/√h | 角度随机游走 |
| 刻度因数K | 1.000372 | 实测值 |
| 零偏稳定性 | ~0.01 deg/h | 短期 |
| 采样率 | 400 Hz | |

### 3.2 原子陀螺 (仿真)
| 参数 | 值 | 说明 |
|------|-----|------|
| ARW | 2e-4 deg/√h | 比FOG好15倍 |
| T_cycle | 2 s | 测量周期 |
| Bias | 1e-5 deg/h | 极小 |

### 3.3 精度对比
```
FOG 1小时角度噪声:  0.003 * √1 = 0.003 deg = 10.8 arcsec
原子 1小时角度噪声: 0.0002 * √1 = 0.0002 deg = 0.72 arcsec
                                            比FOG好15倍
```

## 4. 实验结果

### 4.1 三种场景对比 (smooth_demo.cpp)

| 场景 | eb设置 | 漂移率 | 说明 |
|------|--------|--------|------|
| 1 | eb = eb_true (已知真值) | ~0.7 nm/h | 理论极限 |
| 2 | eb = 0 (不补偿) | ~1.5 nm/h | 最差情况 |
| 3 | 回溯平滑 | ~0.8 nm/h | 接近极限 |

### 4.2 关键发现: FOG零偏漂移

分析1小时数据发现，FOG Z轴角速率持续漂移：
```
时间 0min:  FOG_z = 7.959 deg/h
时间 60min: FOG_z = 7.8 deg/h
漂移: -0.16 deg/h (约2.8 arcsec/min)
```

而原子陀螺保持稳定在 ~7.97 deg/h，证明原子陀螺可提供长期稳定基准。

## 5. 策略对比 (strategy_compare.cpp)

### 5.1 测试的6种策略

| 策略 | 方法 | 结果 |
|------|------|------|
| None | 不修正 | 基线 |
| First Window | 用第1个窗口的eb | 初期好，后期差 |
| Per Window | 每窗口独立eb | 最优 |
| Global Mean | 所有窗口均值 | 中等 |
| Linear Fit | 线性拟合eb(t) | 较好 |
| Sliding Average | 滑动平均 | 较好 |

### 5.2 最佳策略: Per Window
每个 T_cycle 独立估计 eb，能跟踪零偏的缓慢漂移。

## 6. 理论分析

### 6.1 原子陀螺角度噪声计算
```
ARW = 2e-4 deg/√h
T_cycle = 2s = 2/3600 h

角度噪声 σ_θ = ARW * √T_cycle
            = 2e-4 * √(2/3600)
            = 2e-4 * 0.0236
            = 4.7e-6 deg
            = 0.017 arcsec (每次测量)
```

### 6.2 eb估计精度
```
在 T_cycle 内:
  FOG积分角度:  θ_fog = (wie_b + eb) * T_cycle + noise_fog
  原子测量角度: θ_cai = wie_b * T_cycle + noise_cai

  eb = (θ_fog - θ_cai) / T_cycle

噪声传播:
  σ_eb = √(σ_fog² + σ_cai²) / T_cycle
       ≈ σ_fog / T_cycle  (原子噪声可忽略)
       = (ARW_fog * √T_cycle) / T_cycle
       = ARW_fog / √T_cycle
       = 0.003 / √(2/3600)
       = 0.127 deg/h (单次估计)

N次平均后:
  σ_eb_avg = 0.127 / √N deg/h
```

### 6.3 提升极限
平滑后的主要误差来源：
1. 原子陀螺测量噪声 (很小)
2. FOG ARW 在每个周期内的累积
3. 姿态误差导致的 wie_b 计算误差

## 7. 代码实现

### 7.1 统一仿真器接口 (cai_sim.h)
```cpp
namespace cai {

struct CAIParams {
    double T_cycle = 2.0;         // s
    double arw_dpsh = 2e-4;       // deg/√h
    double bias_dph = 1e-5;       // deg/h
};

class AtomicGyroSimulator {
public:
    void Init(const Vector3d& att_rad);
    Vector3d Measure(double dt);      // 返回角增量
    Vector3d GetTrueWieB() const;     // 真实 wie_b
    Vector3d GetAngleNoise();         // 噪声
};

}  // namespace cai
```

### 7.2 使用示例
```cpp
#include "cai_sim.h"

// 初始化
cai::AtomicGyroSimulator atom(pos, glv);
atom.Init(att_true);

// 每个周期
Vector3d wie_b_ref = atom.GetTrueWieB();
Vector3d eb_est = (fog_sum / T_cycle) - wie_b_ref;
```

## 8. 重要教训

### 8.1 wie_b 的计算必须准确
```cpp
// 正确: 从真实姿态计算
Matrix3d Cnb_true = ...;  // 已知或高精度估计
Vector3d wie_b = Cnb_true.transpose() * wie_n;

// 错误: 从粗对准姿态计算
Matrix3d Cnb_coarse = ...;  // 有误差
Vector3d wie_b = Cnb_coarse.transpose() * wie_n;  // 会传递误差!
```

### 8.2 回溯 vs 实时
- **回溯平滑**: 可以用"未来"数据，精度高
- **实时滤波**: 只能用当前和过去数据，精度受限

### 8.3 仿真 vs 实测
本实验使用仿真原子陀螺，假设 wie_b 已知。实际应用中：
- 原子陀螺测量本身有不确定性
- 需要考虑启动时间、死区等

## 9. 项目文件结构

### 9.1 核心代码 (src/)
| 文件 | 说明 |
|------|------|
| `cai_sim.h` | 原子陀螺仿真器头文件 |
| `cai_sim.cpp` | 原子陀螺仿真器实现 |
| `eskf_cai.h` | ESKF 15维滤波器 (含CAIGSimulator) |
| `support.h` | 数据加载、辅助函数 |
| `test.cpp` | 主测试程序 |
| `hybirdmain.cpp` | 混合导航主程序 (开发中) |

### 9.2 PSINS库 (psins/)
| 文件 | 说明 |
|------|------|
| `glv.h` | 全局常量 (地球参数、单位转换) |
| `earth.h/cpp` | 地球模型 (重力、曲率半径) |
| `ins_math.h/cpp` | 数学工具 (四元数、旋转矩阵) |
| `ins_state.h/cpp` | INS状态定义 |
| `kf_state.h/cpp` | 卡尔曼滤波器 |
| `sins_engine.h/cpp` | 惯导解算引擎 |

### 9.3 归档代码 (archive/)
| 文件 | 说明 |
|------|------|
| `smooth_demo.cpp` | 回溯平滑演示 (使用cai_sim.h) |
| `strategy_compare.cpp` | 6种平滑策略对比 |
| `fgo_test.cpp` | ⭐ FGO成功案例 (0.72 nm/h) |
| `fgo_sim_test.cpp` | FGO仿真验证 |
| `theory_test.cpp` | 理论极限分析 |
| `example_minimal.cpp` | 最简使用示例 |
| `example_complete.cpp` | 完整使用示例 |

### 9.4 数据文件
| 文件 | 说明 |
|------|------|
| `fog3h.csv` | 3小时FOG静态数据 (400Hz) |
| `fog.csv` | 1小时FOG数据 |

### 9.5 文档 (doc/)
| 文件 | 说明 |
|------|------|
| `optimal_params.md` | 最优对准参数记录 |
| `fgo_experiments.md` | FGO实验总结 |
| `atom_gyro_smoothing.md` | 本文档 |

## 10. 未来工作

1. **实时ESKF**: 将原子陀螺作为观测，实时估计零偏
2. **组合导航**: FOG高频 + 原子低频融合
3. **实测验证**: 用真实原子陀螺数据测试

## 11. 关键结论

1. ✅ 原子陀螺可显著提升零偏估计精度
2. ✅ 回溯平滑是有效的离线处理方法
3. ✅ Per-window 策略最适合跟踪漂移
4. ⚠️ 核心难点是准确的 wie_b 计算
5. ❌ 直接在FGO中添加原子约束容易失败
