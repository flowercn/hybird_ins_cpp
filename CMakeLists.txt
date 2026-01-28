cmake_minimum_required(VERSION 3.10)
project(HybridINS)

# --- 0. 基础设置 ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release) # 建议使用 Release 以获得最快计算速度

# --- 1. 寻找依赖库 ---
find_package(Eigen3 REQUIRED)
find_package(GTSAM REQUIRED)
find_package(Boost REQUIRED)

# --- 2. 全局头文件路径配置 ---
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/psins
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${EIGEN3_INCLUDE_DIR}
    ${GTSAM_INCLUDE_DIR}
)

# --- 3. 定义源文件与库 ---

# [核心库] PSINS (将底层算法封装为静态库)
set(PSINS_SRC 
    psins/earth.cpp 
    psins/ins_math.cpp 
    psins/ins_state.cpp 
    psins/kf_state.cpp
    psins/sins_engine.cpp
)
add_library(psins STATIC ${PSINS_SRC})
target_link_libraries(psins Eigen3::Eigen)

# FGO 模块源文件
set(FGO_SRC
    src/fgo_core.cpp
    src/cai_sim.cpp
)

# --- 4. 编译目标 ---

# [目标 A] 主程序 fgo_main
add_executable(fgo_main 
    src/fgo_main.cpp 
    ${FGO_SRC}
)
target_link_libraries(fgo_main 
    psins
    gtsam 
    gtsam_unstable
    ${Boost_LIBRARIES}
)

# [目标 B] 调试程序 test
# 修正：移除不存在的 src/support.cpp，只包含头文件供 IDE 索引（可选）
add_executable(test src/test.cpp src/support.h)
target_link_libraries(test psins)

# [目标 C] 平滑演示程序 smooth_demo
# 修正：移除不存在的 src/support.cpp
add_executable(smooth_demo src/smooth_demo.cpp src/support.h)
target_link_libraries(smooth_demo psins)