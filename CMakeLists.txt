cmake_minimum_required(VERSION 3.0.2)
project(test_imu)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=pedantic")

find_package(Eigen3 REQUIRED)
find_package(catkin REQUIRED COMPONENTS
  roscpp
  srrg2_core
  srrg2_core_ros
  srrg2_solver
  srrg_cmake_modules
)

catkin_package(
  INCLUDE_DIRS
    src
  CATKIN_DEPENDS
    srrg_cmake_modules
    srrg2_core
    srrg2_core_ros
    srrg2_solver
  DEPENDS
    EIGEN3
)

include(${srrg_cmake_modules_INCLUDE_DIRS}/CMakeCompileOptions.txt)

include_directories(
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/src
)

# this is important. it tells cmake where to look for libraries
link_directories(${catkin_LIBRARY_DIRS})

add_subdirectory(src)
add_subdirectory(examples)