add_definitions(-DNCLT_DATA_FOLDER="${PROJECT_SOURCE_DIR}/data/nclt")

add_executable(imu_gps_synthetic imu_gps_synthetic.cpp)
target_link_libraries(imu_gps_synthetic
                      common
                      imu_preintegrator
                      synthetic
                      preintegration_factor
                      srrg2_solver_core_library
                      srrg2_solver_types3d_library
)

add_executable(view_graph view_graph.cpp)
target_link_libraries(view_graph
                      common
                      imu_preintegrator
                      preintegration_factor
                      srrg2_solver_core_library
                      srrg2_solver_types3d_library
                      ${QGLViewer_LIBRARIES}
                      ${SRRG_QT_LIBRARIES}
                      ${GLUT_LIBRARIES}
                      ${catkin_LIBRARIES}
                      GL
                      pthread
)

add_executable(synthetic_gtsam synthetic_gtsam.cpp)
target_link_libraries(synthetic_gtsam
                      gtsam
                      synthetic
)

find_package(GTSAMCMakeTools)
find_package(GTSAM REQUIRED) # Uses installed package
include_directories(${GTSAM_INCLUDE_DIR})