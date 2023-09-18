#pragma once

#include <Eigen/Dense>
#include <srrg_geometry/geometry_defs.h>

namespace test_imu {

  namespace core = srrg2_core;

  struct ImuMeasurement {
    core::Vector3f acceleration;
    core::Vector3f angular_vel;

    double timestamp;
  };

} // namespace test_imu