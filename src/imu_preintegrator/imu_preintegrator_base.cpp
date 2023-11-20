#include "imu_preintegrator_base.h"

namespace test_imu {

  ImuPreintegratorBase::ImuPreintegratorBase() {
    reset();
  }

  void ImuPreintegratorBase::reset() {
    measurements_.clear();
    dT_ = 0;

    // nominal values for the bias
    bias_acc_  = Vector3::Zero();
    bias_gyro_ = Vector3::Zero();

    // bias correction
    dR_db_gyro_ = Matrix3::Zero();
    dv_db_acc_  = Matrix3::Zero();
    dv_db_gyro_ = Matrix3::Zero();
    dp_db_acc_  = Matrix3::Zero();
    dp_db_gyro_ = Matrix3::Zero();
  }

  const void ImuPreintegratorBase::setNoiseGyroscope(const core::Vector3f& v) {
    sigma_noise_.block<3, 3>(NGidx, NGidx) = v.cast<Scalar>().asDiagonal();
  }
  const void ImuPreintegratorBase::setNoiseAccelerometer(const core::Vector3f& v) {
    sigma_noise_.block<3, 3>(NAidx, NAidx) = v.cast<Scalar>().asDiagonal();
  }
  const void ImuPreintegratorBase::setNoiseBiasGyroscope(const core::Vector3f& v) {
    sigma_noise_.block<3, 3>(NBGidx, NBGidx) = v.cast<Scalar>().asDiagonal();
  }
  const void ImuPreintegratorBase::setNoiseBiasAccelerometer(const core::Vector3f& v) {
    sigma_noise_.block<3, 3>(NBAidx, NBAidx) = v.cast<Scalar>().asDiagonal();
  }

  void ImuPreintegratorBase::getPrediction(const core::Isometry3f& Ti,
                                           const core::Vector3f& vi,
                                           core::Isometry3f& Tf,
                                           core::Vector3f& vf) {
    Tf.setIdentity();

    float T = measurements_.back().timestamp - measurements_.at(0).timestamp;

    vf = Ti.linear() * delta_v() + vi;

    Tf.linear()      = Ti.linear() * delta_R();
    Tf.translation() = Ti.linear() * delta_p() + Ti.translation() + T * vi;
  }
} // namespace test_imu