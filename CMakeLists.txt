cmake_minimum_required(VERSION 3.10)
project(HybridINS)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

# --- 1. 寻找依赖 ---
find_package(Eigen3 REQUIRED)
find_package(GTSAM REQUIRED)
find_package(Boost REQUIRED)

# --- 2. [关键] 全局头文件路径配置 ---
# 这样配置后:
# 1. 引用 psins 库文件直接用 #include "glv.h"
# 2. 引用 src 下文件直接用 #include "fgo_core.h"
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/psins
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${EIGEN3_INCLUDE_DIR}
    ${GTSAM_INCLUDE_DIR} # 确保 GTSAM 头文件被包含
)

# --- 3. 定义源文件集合 ---
set(PSINS_SRC 
    psins/glv.cpp 
    psins/earth.cpp 
    psins/ins_math.cpp 
    psins/ins_state.cpp 
    psins/kf_state.cpp
    psins/sins_engine.cpp
)

set(FGO_SRC
    src/fgo_core.cpp
    src/cai_sim.cpp
)

# --- 4. 编译主程序 fgo_main ---
add_executable(fgo_main 
    src/fgo_main.cpp 
    ${PSINS_SRC} 
    ${FGO_SRC}
)

target_link_libraries(fgo_main 
    Eigen3::Eigen 
    gtsam 
    ${Boost_LIBRARIES}
)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/test_gtsam_logic.cpp")
    add_executable(test_gtsam src/test_gtsam_logic.cpp)
    target_link_libraries(test_gtsam gtsam)
endif()
