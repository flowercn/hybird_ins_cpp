import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import re

# ==========================================
# [配置] 绘图风格
# ==========================================
def set_paper_style():
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['axes.titlesize'] = 16
    plt.rcParams['legend.fontsize'] = 11
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    plt.rcParams['lines.linewidth'] = 2.0
    plt.rcParams['axes.grid'] = True
    plt.rcParams['grid.alpha'] = 0.4
    plt.rcParams['grid.linestyle'] = '--'
    plt.rcParams['axes.unicode_minus'] = False 
    
    fonts = ['WenQuanYi Micro Hei', 'Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    for f in fonts:
        try:
            plt.rcParams['font.sans-serif'] = [f]
            break
        except:
            continue
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['#D62728', '#1F77B4', '#2CA02C', '#FF7F0E', '#9467BD'])

# ==========================================
# [核心] 数据清洗
# ==========================================
def process_dataframe(df, cycle_time=2.0):
    df.columns = df.columns.str.strip()
    
    # 智能降采样
    mask = (np.abs(df['time'] % cycle_time) < 0.01) | (np.abs(df['time'] % cycle_time - cycle_time) < 0.01)
    df_resampled = df[mask].drop_duplicates(subset='time').copy()
    
    # -----------------------------------------------------------
    # [修改点] 1. 陀螺累计残差 (弧度 -> 角秒 ″)
    # -----------------------------------------------------------
    if 'res_gyro_x' in df_resampled.columns:
        # 1 rad = (180/pi) * 3600 arcsec
        rad2arcsec = (180.0 / np.pi) * 3600.0
        df_resampled['cum_gyro_x'] = df_resampled['res_gyro_x'].cumsum() * rad2arcsec
        df_resampled['cum_gyro_y'] = df_resampled['res_gyro_y'].cumsum() * rad2arcsec
        df_resampled['cum_gyro_z'] = df_resampled['res_gyro_z'].cumsum() * rad2arcsec
    
    # 2. 加计累计残差 (m/s)
    if 'res_acc_x' in df_resampled.columns:
        df_resampled['cum_acc_x'] = df_resampled['res_acc_x'].cumsum()
        df_resampled['cum_acc_y'] = df_resampled['res_acc_y'].cumsum()
        df_resampled['cum_acc_z'] = df_resampled['res_acc_z'].cumsum()

    return df, df_resampled

# ==========================================
# [主程序] 绘图逻辑
# ==========================================
def main():
    set_paper_style()
    
    search_pattern = "nav_dead_*.csv"
    files = glob.glob(search_pattern)
    files.sort()
    
    print(f"=== 死区效应分析绘图工具 ===")
    
    if not files:
        print(f"未找到数据文件 ({search_pattern})")
        return

    datasets = []
    print(f"正在处理 {len(files)} 个文件...")

    for f in files:
        try:
            raw_df = pd.read_csv(f)
            match = re.search(r"nav_dead_(\d+\.?\d*)s", f)
            if match:
                t_active = float(match.group(1))
                t_dead = 2.0 - t_active
                label = f"死区 {t_dead:.1f}s (工作 {t_active}s)"
                sort_key = t_dead
            else:
                label = f
                sort_key = 0
            
            full_df, resampled_df = process_dataframe(raw_df)
            
            datasets.append({
                'label': label,
                'full': full_df,
                'resampled': resampled_df,
                'sort': sort_key
            })
            
        except Exception as e:
            print(f"处理 {f} 失败: {e}")

    datasets.sort(key=lambda x: x['sort'], reverse=True)
    
    output_dir = "figs_dead_time"
    if not os.path.exists(output_dir): os.makedirs(output_dir)

    # 图 1: 位置误差对比 (不变)
    fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for ds in datasets:
        df = ds['full']
        t = df['time'] / 3600.0
        ax1.plot(t, df['lat_err'], label=ds['label'])
        ax2.plot(t, df['lon_err'], label=ds['label'])
    
    ax1.set_title('位置误差对比')
    ax1.set_ylabel('北向误差 (m)')
    ax1.legend(loc='upper left', ncol=2)
    ax2.set_ylabel('东向误差 (m)')
    ax2.set_xlabel('时间 (h)')
    ax2.set_xlim(0, 24)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/1_位置误差.png", dpi=300)
    print(f"已生成: {output_dir}/1_位置误差.png")

    # -----------------------------------------------------------
    # [修改点] 图 2: 陀螺累计残差 (单位改成 ″)
    # -----------------------------------------------------------
    if 'cum_gyro_x' in datasets[0]['resampled'].columns:
        fig2, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
        
        for ds in datasets:
            df = ds['resampled']
            t = df['time'] / 3600.0
            axes[0].plot(t, df['cum_gyro_x'], label=ds['label'])
            axes[1].plot(t, df['cum_gyro_y'], label=ds['label'])
            axes[2].plot(t, df['cum_gyro_z'], label=ds['label'])
            
        axes[0].set_title('陀螺累计残差')
        axes[0].set_ylabel('X轴角度漂移 (″)')  # 改为角秒
        axes[0].legend(loc='upper left')
        
        axes[1].set_ylabel('Y轴角度漂移 (″)')  # 改为角秒
        axes[2].set_ylabel('Z轴角度漂移 (″)')  # 改为角秒
        axes[2].set_xlabel('时间 (h)')
        axes[2].set_xlim(0, 24)
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/2_陀螺残差.png", dpi=300)
        print(f"已生成: {output_dir}/2_陀螺残差.png")

    # 图 3: 加计累计残差 (不变)
    if 'cum_acc_x' in datasets[0]['resampled'].columns:
        fig3, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
        for ds in datasets:
            df = ds['resampled']
            t = df['time'] / 3600.0
            axes[0].plot(t, df['cum_acc_x'], label=ds['label'])
            axes[1].plot(t, df['cum_acc_y'], label=ds['label'])
            axes[2].plot(t, df['cum_acc_z'], label=ds['label'])
            
        axes[0].set_title('加计累计残差')
        axes[0].set_ylabel('X轴速度漂移 (m/s)')
        axes[0].legend(loc='upper left')
        axes[1].set_ylabel('Y轴速度漂移 (m/s)')
        axes[2].set_ylabel('Z轴速度漂移 (m/s)')
        axes[2].set_xlabel('时间 (h)')
        axes[2].set_xlim(0, 24)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/3_加计残差.png", dpi=300)
        print(f"已生成: {output_dir}/3_加计残差.png")

if __name__ == "__main__":
    main()