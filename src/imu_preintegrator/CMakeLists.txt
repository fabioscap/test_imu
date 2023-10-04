add_library(imu_preintegrator SHARED
imu_preintegrator_slim.cpp
imu_preintegrator.cpp)

target_link_libraries(imu_preintegrator
                      common)