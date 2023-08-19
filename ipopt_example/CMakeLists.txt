project(ipopt_example)

cmake_minimum_required (VERSION 3.5)

# -g allows for gdb debugging
# turn on -03 for best performance
add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS "${CXX_FLAGS}")

set(sources src/example.cpp src/get_started.cpp)

include_directories(
    include
)

add_executable(cppad_ipopt ${sources})
target_link_libraries(cppad_ipopt ipopt)