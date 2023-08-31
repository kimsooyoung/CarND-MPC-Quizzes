project(mpc)

cmake_minimum_required (VERSION 3.5)

# -g allows for gdb debugging
# turn on -03 for best performance
add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS "${CXX_FLAGS}")

# Find the Python3 package
find_package(Eigen3 REQUIRED)
find_package(Python3 REQUIRED)
find_package (Python3 COMPONENTS Interpreter Development)

# include_directories(src/Eigen-3.3)
include_directories(
    include
    ${EIGEN3_INCLUDE_DIRS}
)

add_executable(circular_traj src/circular_traj.cpp)
target_include_directories(circular_traj PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(circular_traj PRIVATE Python3::Python)

add_executable(mpc_traj src/main.cpp src/MPC.cpp)
target_include_directories(mpc_traj PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(mpc_traj PRIVATE Python3::Python ipopt)
