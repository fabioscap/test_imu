add_library(preintegration_factor SHARED
            imu_preintegration_factor.cpp
            )

target_link_libraries(preintegration_factor PUBLIC
  srrg2_solver_core_library
  srrg2_solver_types3d_library
  common
  ${catkin_LIBRARIES}
)