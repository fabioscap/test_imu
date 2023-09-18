add_library(variables_and_factors
            imu_preintegration_factor.h
            imu_preintegration_factor.cpp)

target_link_libraries(variables_and_factors
                      srrg2_solver_core_library
                      common)