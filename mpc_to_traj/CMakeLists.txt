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

add_executable(line_traj src/line_traj.cpp)
target_include_directories(line_traj PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(line_traj PRIVATE Python3::Python)

add_executable(l_shape_traj src/l_shape_traj.cpp)
target_include_directories(l_shape_traj PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(l_shape_traj PRIVATE Python3::Python)

add_executable(square_shape_traj src/square_shape_traj.cpp)
target_include_directories(square_shape_traj PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(square_shape_traj PRIVATE Python3::Python)

add_executable(mpc_traj_line src/line_main.cpp src/MPC.cpp)
target_include_directories(mpc_traj_line PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(mpc_traj_line PRIVATE Python3::Python ipopt)

add_executable(mpc_traj_line_qd src/line_main_qd.cpp src/MPCQD.cpp)
target_include_directories(mpc_traj_line_qd PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(mpc_traj_line_qd PRIVATE Python3::Python ipopt)

add_executable(mpc_traj_circle src/circle_main.cpp src/MPC.cpp)
target_include_directories(mpc_traj_circle PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(mpc_traj_circle PRIVATE Python3::Python ipopt)

add_executable(square_main_qd src/square_main_qd.cpp src/MPCQD.cpp)
target_include_directories(square_main_qd PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(square_main_qd PRIVATE Python3::Python ipopt)
