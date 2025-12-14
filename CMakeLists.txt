cmake_minimum_required(VERSION 3.10)

# 设置项目名称
project(CircuitSimulator LANGUAGES CXX)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# 默认使用 Release（O3）构建，除非用户显式指定
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(Eigen3 3.4 REQUIRED NO_MODULE)


# 添加源代码文件
file(GLOB SRC_FILES
    src/*.cpp
)
# 添加可执行文件
add_executable(${PROJECT_NAME} ${SRC_FILES})

# Release 下关掉 assert + Eigen debug checks（对 Eigen 很重要）
target_compile_definitions(CircuitSimulator PRIVATE
  $<$<CONFIG:Release>:NDEBUG>
  $<$<CONFIG:Release>:EIGEN_NO_DEBUG>
)

# Release 编译优化：O3 + 本机指令集
target_compile_options(CircuitSimulator PRIVATE
  $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU>>:-O3 -march=native -mtune=native -pipe>
  $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:Clang>>:-O3 -march=native -mtune=native>
)

# LTO/IPO（链接时优化，通常能再快一些）
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_ok OUTPUT ipo_err)
if(ipo_ok)
  set_property(TARGET CircuitSimulator PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()

# 可选：fast-math（可能影响数值稳定/收敛，默认关闭以保持精度）
option(CSIM_FAST_MATH "Enable fast-math (may hurt numerical stability)" OFF)
if(CSIM_FAST_MATH)
  target_compile_options(CircuitSimulator PRIVATE
    $<$<AND:$<CONFIG:Release>,$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>>:-ffast-math -fno-math-errno>
  )
endif()

# 将 include 文件夹加入头文件搜索路径
target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(${PROJECT_NAME}
    PRIVATE
        Eigen3::Eigen
)



target_compile_options(${PROJECT_NAME}
    PRIVATE
        -Wall -Wextra -Wpedantic
)
