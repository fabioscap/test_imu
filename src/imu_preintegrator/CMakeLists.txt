add_library(imu_preintegrator SHARED
imu_preintegrator_base.cpp
imu_preintegrator.cpp
imu_preintegrator_ukf.cpp
)

target_link_libraries(imu_preintegrator PUBLIC
                      common)