#pragma once

#include <Eigen/Dense>
#include <srrg_geometry/geometry_defs.h>

namespace test_imu {

  const srrg2_core::Vector3f gravity_vector(0, 0, -9.81);

  struct IMUMeasurement {
    srrg2_core::Vector3f acceleration;
    srrg2_core::Vector3f angular_vel;

    double timestamp;
  };

  class IMUSensor {
  public:
    IMUSensor(float freq);

  protected:
    float freq_; // the frequency at which a new measurement becomes available
  };

} // namespace test_imu