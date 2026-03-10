cmake_minimum_required(VERSION 3.10)
project(HybridINS)

# --- 0. 基础设置 ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# 核心实验跑数据，必须使用 Release 以获得最高执行效率
set(CMAKE_BUILD_TYPE Release) 

# --- 1. 寻找依赖库 ---
# 既然是纯解析推导，我们只需要极其轻量的 Eigen3 进行矩阵运算
find_package(Eigen3 REQUIRED)

# --- 2. 全局头文件路径配置 ---
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/psins
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${EIGEN3_INCLUDE_DIR}
)

# --- 3. 定义源文件与核心库 ---

# [核心库] PSINS (惯导底层引擎)
set(PSINS_SRC 
    psins/earth.cpp 
    psins/ins_math.cpp 
    psins/ins_state.cpp 
    psins/kf_state.cpp
    psins/sins_engine.cpp
)
add_library(psins STATIC ${PSINS_SRC})
target_link_libraries(psins Eigen3::Eigen)

# --- 4. 编译目标 ---

# [唯一主程序] 顶刊论文核心算法对决 (test_fgo)
add_executable(test_fgo src/test_fgo.cpp src/experiment_app.cpp src/cai_sim.cpp)
target_link_libraries(test_fgo psins)
