import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# --- 常数定义 (与 GLV 保持一致) ---
Re = 6378137.0
deg = np.pi / 180.0

# --- 绘图风格设置 ---
plt.style.use('seaborn-v0_8-paper')
# plt.rcParams['font.sans-serif'] = ['SimHei'] # 用来正常显示中文标签 (如果系统不支持可注释掉)
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

def xygo(ax, ylabel, xlabel='Time (s)'):
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.legend(loc='best')

def plot_results():
    # 1. 文件读取 (兼容相对路径和绝对路径)
    file_path = 'res_sins_engine.csv'
    
    # 如果当前目录没找到，尝试去 build 目录找 (方便在项目根目录运行)
    if not os.path.exists(file_path):
        file_path = '../build/res_sins_engine.csv'
    
    # 再次检查
    if not os.path.exists(file_path):
        # 最后尝试用户提供的绝对路径作为 fallback
        file_path = '/home/v/dev/hybrid_ins_cpp/build/res_sins_engine.csv'
    
    try:
        data = pd.read_csv(file_path)
        print(f"Successfully loaded data from: {file_path}")
        print(f"Data columns: {data.columns.tolist()}")
    except FileNotFoundError:
        print(f"Error: Could not find 'res_sins_engine.csv' in current, ../build/, or absolute path.")
        return

    t = data['t'].values
    
    # 初始位置 (参考 C++ 代码中的 cfg.pos_ref)
    # 注意：这里我们用数据的第一行作为基准，或者手动指定
    lat0 = data['lat'].iloc[0] # 32.0286
    lon0 = data['lon'].iloc[0] # 118.8533
    h0   = data['h'].iloc[0]   # 17.0

    # ==========================================
    # 图 1: 导航状态 (速度 + 位置误差)
    # ==========================================
    fig1 = plt.figure(figsize=(12, 10))
    
    # --- 1.1 速度 ---
    ax1 = fig1.add_subplot(3, 1, 1)
    ax1.plot(t, data['vE'], label='vE')
    ax1.plot(t, data['vN'], label='vN')
    ax1.plot(t, data['vU'], label='vU (Should be 0)', linewidth=2, linestyle='--')
    ax1.set_title('Velocity (m/s) [Includes Vertical Damping Check]')
    xygo(ax1, 'Vel (m/s)')

    # --- 1.2 水平位置漂移 (转换为米) ---
    # dLat = (Lat - Lat0) * Re
    # dLon = (Lon - Lon0) * Re * cos(Lat0)
    d_lat_m = (data['lat'] - lat0) * deg * Re
    d_lon_m = (data['lon'] - lon0) * deg * Re * np.cos(lat0 * deg)
    
    ax2 = fig1.add_subplot(3, 1, 2)
    ax2.plot(t, d_lon_m, label='East Drift (m)')
    ax2.plot(t, d_lat_m, label='North Drift (m)')
    ax2.set_title('Horizontal Position Drift')
    xygo(ax2, 'Drift (m)')

    # --- 1.3 高度漂移 ---
    d_h_m = data['h'] - h0
    ax3 = fig1.add_subplot(3, 1, 3)
    ax3.plot(t, d_h_m, label='Height Drift', color='green')
    ax3.set_title('Vertical Position Drift (Should be locked)')
    xygo(ax3, 'dH (m)')

    plt.tight_layout()
    plt.savefig('sins_nav_status.png')

    # ==========================================
    # 图 2: 姿态角 (Attitude)
    # ==========================================
    fig2 = plt.figure(figsize=(10, 6))
    ax4 = fig2.add_subplot(1, 1, 1)
    ax4.plot(t, data['pitch'], label='Pitch')
    ax4.plot(t, data['roll'], label='Roll')
    ax4.plot(t, data['yaw'], label='Yaw')
    ax4.set_title('Attitude (deg)')
    xygo(ax4, 'Angle (deg)')
    
    plt.tight_layout()
    plt.savefig('sins_attitude.png')

    # ==========================================
    # 图 3: 零偏曲线 (Bias) [新增]
    # ==========================================
    # 检查数据中是否存在 ebX 列
    if 'ebX' in data.columns:
        fig3 = plt.figure(figsize=(12, 8))
        
        # 陀螺零偏
        ax5 = fig3.add_subplot(2, 1, 1)
        ax5.plot(t, data['ebX'], label='ebX')
        ax5.plot(t, data['ebY'], label='ebY')
        ax5.plot(t, data['ebZ'], label='ebZ')
        ax5.set_title('Gyro Bias Estimate (deg/h)')
        xygo(ax5, 'Bias (deg/h)')

        # 加计零偏
        ax6 = fig3.add_subplot(2, 1, 2)
        ax6.plot(t, data['dbX'], label='dbX')
        ax6.plot(t, data['dbY'], label='dbY')
        ax6.plot(t, data['dbZ'], label='dbZ')
        ax6.set_title('Accel Bias Estimate (ug)')
        xygo(ax6, 'Bias (ug)')

        plt.tight_layout()
        plt.savefig('sins_biases.png')
        print("Saved bias plot to sins_biases.png")

    # ==========================================
    # 分段验证检查 (Separated Validation)
    # ==========================================
    align_time = 300.0
    
    # 1. 对准阶段数据 (t <= 300)
    mask_align = t <= align_time
    max_vu_align = np.max(np.abs(data['vU'][mask_align])) if np.any(mask_align) else 0.0
    
    # 2. 导航阶段数据 (t > 300)
    mask_nav = t > align_time
    if np.any(mask_nav):
        # 导航阶段的垂直速度（理论上应为0）
        max_vu_nav = np.max(np.abs(data['vU'][mask_nav]))
        # 导航阶段的高度漂移（相对于 t=300 时刻的高度变化，或者相对于 h0）
        max_dh_nav = np.max(np.abs(d_h_m[mask_nav]))
    else:
        max_vu_nav = 0.0
        max_dh_nav = 0.0

    print("-" * 50)
    print(f"📊 [SINS Phase Analysis Report]")
    print(f"   Total Duration:    {t[-1]:.2f} s")
    print(f"   Alignment Time:    {align_time} s")
    print("-" * 50)
    print(f"🛑 [Alignment Phase] (KF Estimation Noise)")
    print(f"   Max Vertical Vel:  {max_vu_align:.6f} m/s (Normal fluctuation)")
    print("-" * 50)
    print(f"🚀 [Navigation Phase] (Pure INS + Damping)")
    print(f"   Max Vertical Vel:  {max_vu_nav:.9f} m/s (Target: 0.0)")
    print(f"   Max Height Drift:  {max_dh_nav:.9f} m   (Target: 0.0)")
    print("-" * 50)
    print(f"📍 [Final State]")
    print(f"   Pos: Lat={data['lat'].iloc[-1]:.6f}, Lon={data['lon'].iloc[-1]:.6f}")
    if 'ebX' in data.columns:
        print(f"   Est Gyro Bias: [{data['ebX'].iloc[-1]:.3f}, {data['ebY'].iloc[-1]:.3f}, {data['ebZ'].iloc[-1]:.3f}] deg/h")
        print(f"   Est Acc Bias:  [{data['dbX'].iloc[-1]:.1f}, {data['dbY'].iloc[-1]:.1f}, {data['dbZ'].iloc[-1]:.1f}] ug")
    print("-" * 50)
    
    plt.show()

if __name__ == "__main__":
    plot_results()