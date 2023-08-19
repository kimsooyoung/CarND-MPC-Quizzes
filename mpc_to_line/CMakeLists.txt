project(mpc)

cmake_minimum_required (VERSION 3.5)

# -g allows for gdb debugging
# turn on -03 for best performance
add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS "${CXX_FLAGS}")

set(sources src/MPC.cpp src/main.cpp)

# Find the Python3 package
find_package(Eigen3 REQUIRED)
find_package(Python3 REQUIRED)
find_package (Python3 COMPONENTS Interpreter Development)

# include_directories(src/Eigen-3.3)
include_directories(
    include
    ${EIGEN3_INCLUDE_DIRS}
)

add_executable(mpc ${sources})
target_include_directories(mpc PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(mpc PRIVATE Python3::Python ipopt)

add_executable(mpc_traj src/MPC_traj.cpp src/main_traj.cpp)
target_include_directories(mpc_traj PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(mpc_traj PRIVATE Python3::Python ipopt)

add_executable(sys_info src/sys_info.cpp)
