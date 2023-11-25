add_executable(integrate_trajectory integrate_trajectory.cpp)
target_link_libraries(integrate_trajectory 
                      srrg2_solver_core_library
                      synthetic
                      common
                      imu_preintegrator
)

add_executable(test_ test.cpp)
target_link_libraries(test_ 
                      srrg2_solver_core_library
                      synthetic
                      common
                      imu_preintegrator
)

add_executable(test_single_factor test_single_factor.cpp)
target_link_libraries(test_single_factor
                      common
                      synthetic
                      imu_preintegrator
                      preintegration_factor
                      srrg2_solver_core_library
                      srrg2_solver_types3d_library
                      ${catkin_LIBRARIES}
)

add_executable(example_kitti_imu_gps example_kitti_imu_gps.cpp)
target_link_libraries(example_kitti_imu_gps
                      common
                      imu_preintegrator
                      preintegration_factor
                      srrg2_solver_core_library
                      srrg2_solver_types3d_library
)
add_executable(imu_gps_synthetic imu_gps_synthetic.cpp)
target_link_libraries(imu_gps_synthetic
                      common
                      imu_preintegrator
                      synthetic
                      preintegration_factor
                      srrg2_solver_core_library
                      srrg2_solver_types3d_library
)

add_executable(example_kitti_imu_gps_gtsam example_kitti_imu_gps_gtsam.cpp)
target_link_libraries(example_kitti_imu_gps_gtsam
                      gtsam
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


add_executable(test_ukf test_ukf.cpp)
target_link_libraries(test_ukf
                      common
)

add_executable(nclt_gtsam nclt_gtsam.cpp)
target_link_libraries(nclt_gtsam
                      gtsam
)

add_executable(nclt nclt.cpp)
target_link_libraries(nclt
                      common
                      imu_preintegrator
                      preintegration_factor
                      srrg2_solver_core_library
                      srrg2_solver_types3d_library
                      ${catkin_LIBRARIES}
)



find_package(GTSAMCMakeTools)
find_package(GTSAM REQUIRED) # Uses installed package
include_directories(${GTSAM_INCLUDE_DIR})