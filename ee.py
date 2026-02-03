import numpy as np

input_file = "/home/v/dev/hybrid_ins_cpp/fog30h.csv"

fs = 400  # 采样率
duration_sec = 6 * 60 * 60  # 6小时
num_rows = fs * duration_sec  # 8640000

data = np.loadtxt(input_file, delimiter=',', max_rows=num_rows)
vec = data[:, 3:6].mean(axis=0)
norm = np.linalg.norm(vec)
unit_vec = vec / norm

np.set_printoptions(precision=12, suppress=False)
print("前三列均值向量:", np.array2string(vec, separator=', '))
print("单位向量:", np.array2string(unit_vec, separator=', '))