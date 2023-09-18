#include "imu_preintegration_factor.h"

#include <srrg_geometry/geometry3d.h>

using namespace test_imu;

void ImuPreintegrator::preintegrate(const ImuMeasurement& m_new) {
  // perform preintegration up until m.timestamp
  // therefore the measurement m is used in the next
  // call to preintegrate (maybe we should use the previous dt?)
  if (measurements_.size() == 0) {
    throw std::runtime_error("ImuPreintegrator::preintegrate| call reset once before preintegrate");
  }
  ImuMeasurement& m = measurements_.back();
  measurements_.push_back(m_new);

  float dt = m.timestamp - t_;
  if (dt < 0)
    throw std::runtime_error("vaffanculo");

  // correct the measurements
  core::Vector3f acc_c     = m.acceleration - bias_acc_;
  core::Vector3f ang_vel_c = m.angular_vel - bias_gyro_;

  // update the position
  auto dv = delta_R_ * (acc_c) *dt;
  delta_p_ += (delta_v_ + 0.5 * dv) * dt;

  // update the velocity
  delta_v_ += dv;

  auto Jr = core::geometry3d::expMapSO3(static_cast<core::Vector3f>(dt * ang_vel_c));

  // update orientation
  delta_R_ *= Jr;

  core::fixRotation(delta_R_);

  // covariance matrix update through iterative first order propagation
  A_.block<3, 3>(0, 0) = delta_R_.transpose();

  auto tmp = -delta_R_ * dt * core::geometry3d::skew(acc_c);

  // A_ ...
  A_.block<3, 3>(3, 0) = tmp;
  A_.block<3, 3>(6, 0) = tmp * 0.5 * dt;
  A_.block<3, 3>(6, 3) *= dt;

  // B_ ...
  B_.block<3, 3>(0, 0) = Jr * dt;
  B_.block<3, 3>(3, 3) = delta_R_ * dt;
  B_.block<3, 3>(6, 3) = 0.5 * delta_R_ * dt * dt;

  sigma_.block<9, 9>(0, 0).noalias() =
    A_ * sigma_.block<9, 9>(0, 0) * A_.transpose() + B_ * sigma_imu_ * B_.transpose();
}

void test_imu::ImuPreintegrator::reset(const ImuMeasurement& measurement) {
  std::cout << "ano\n";
}
