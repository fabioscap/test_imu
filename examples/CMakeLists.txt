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