add_library(preintegration_factor SHARED
            imu_preintegration_factor.cpp
            instances.cpp
            variable_se3_expmap_right.cpp
            gps_factor_ad.cpp
            )

target_link_libraries(preintegration_factor PUBLIC
  srrg2_solver_core_library
  srrg2_solver_types3d_library
  imu_preintegrator
  ${catkin_LIBRARIES}
)