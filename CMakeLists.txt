cmake_minimum_required(VERSION 3.0.2)
project(test_imu)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=pedantic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

find_package(Eigen3 REQUIRED)
find_package(catkin REQUIRED COMPONENTS
  roscpp
  srrg2_core
  srrg2_core_ros
  srrg2_solver
  srrg_cmake_modules
  srrg2_qgl_viewport
)

include(${srrg_cmake_modules_INCLUDE_DIRS}/CMakeCompileOptions.txt)

find_package(QGLViewer REQUIRED)
find_package(GLUT REQUIRED)

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

include_directories(
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/src
  ${QGLViewer_INCLUDE_DIR}
  ${SRRG_QT_INCLUDE_DIRS}
  ${GLUT_INCLUDE_DIRS}
)

# this is important. it tells cmake where to look for libraries
link_directories(${catkin_LIBRARY_DIRS})


add_definitions(-DPROJECT_FOLDER="${PROJECT_SOURCE_DIR}")


add_subdirectory(src)
add_subdirectory(examples)