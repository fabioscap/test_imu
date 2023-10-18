#pragma once

#include <Eigen/Dense>
#include <srrg_geometry/geometry_defs.h>

#include <memory>

namespace test_imu {

  // single or double precision for delta increment computation and covariance propagation
  using Scalar = double;

  namespace core = srrg2_core;

  struct ImuMeasurement {
    core::Vector3f acceleration;
    core::Vector3f angular_vel;

    double timestamp;
  };

  template <typename Scalar>
  Eigen::Matrix<Scalar, 3, 3> Rx(Scalar angle_radians) {
    Scalar cos_angle = std::cos(angle_radians);
    Scalar sin_angle = std::sin(angle_radians);

    Eigen::Matrix<Scalar, 3, 3> rotation_matrix;
    rotation_matrix << 1, 0, 0, 0, cos_angle, -sin_angle, 0, sin_angle, cos_angle;
    return rotation_matrix;
  }

  template <typename Scalar>
  Eigen::Matrix<Scalar, 3, 3> Ry(Scalar angle_radians) {
    Scalar cos_angle = std::cos(angle_radians);
    Scalar sin_angle = std::sin(angle_radians);

    Eigen::Matrix<Scalar, 3, 3> rotation_matrix;
    rotation_matrix << cos_angle, 0, sin_angle, 0, 1, 0, -sin_angle, 0, cos_angle;
    return rotation_matrix;
  }

  template <typename Scalar>
  Eigen::Matrix<Scalar, 3, 3> Rz(Scalar angle_radians) {
    Scalar cos_angle = std::cos(angle_radians);
    Scalar sin_angle = std::sin(angle_radians);

    Eigen::Matrix<Scalar, 3, 3> rotation_matrix;
    rotation_matrix << cos_angle, -sin_angle, 0, sin_angle, cos_angle, 0, 0, 0, 1;
    return rotation_matrix;
  }

} // namespace test_imu