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


find_package(GTSAMCMakeTools)
find_package(GTSAM REQUIRED) # Uses installed package
include_directories(${GTSAM_INCLUDE_DIR})