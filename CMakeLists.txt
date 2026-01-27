cmake_minimum_required(VERSION 3.10)
project(HybridINS)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)

# --- 1. 寻找依赖库 ---
find_package(Eigen3 REQUIRED)
find_package(GTSAM REQUIRED)
find_package(Boost REQUIRED)

# --- 2. 全局头文件路径配置 ---
# 这样配置后，源码中 #include "glv.h" 或 #include "support.h" 都能找到
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/psins
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${EIGEN3_INCLUDE_DIR}
    ${GTSAM_INCLUDE_DIR}
)

# --- 3. 定义源文件集合 ---

# [修正] PSINS 模块源文件 (移除了 glv.cpp)
set(PSINS_SRC 
    psins/earth.cpp 
    psins/ins_math.cpp 
    psins/ins_state.cpp 
    psins/kf_state.cpp
    psins/sins_engine.cpp
)

# FGO 模块源文件
set(FGO_SRC
    src/fgo_core.cpp
    src/cai_sim.cpp
)

# --- 4. 编译目标 ---

# [目标 A] 主程序 fgo_main
# 包含 PSINS 和 FGO 所有源码
add_executable(fgo_main 
    src/fgo_main.cpp 
    ${PSINS_SRC} 
    ${FGO_SRC}
)

# 链接 GTSAM 和 Boost
target_link_libraries(fgo_main 
    Eigen3::Eigen 
    gtsam 
    gtsam_unstable
    ${Boost_LIBRARIES}
)

# [目标 B] 调试程序 test
# 目前仅用于调试 PSINS 引擎，因此只编译 PSINS 源码
# 如果 test.cpp 后来引入了 fgo_core 或 cai_sim，需要把 ${FGO_SRC} 加进去
add_executable(test src/test.cpp ${PSINS_SRC})

# test 目前只需要 Eigen，如果后续加入 GTSAM 相关代码，需参考 fgo_main 添加链接库
target_link_libraries(test Eigen3::Eigen)