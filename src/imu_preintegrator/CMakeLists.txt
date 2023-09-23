add_library(imu_preintegrator
imu_preintegrator.h
imu_preintegrator.cpp)

target_link_libraries(imu_preintegrator
                      common)