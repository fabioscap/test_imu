#include "imu_preintegrator_base.h"

namespace test_imu {

  void ImuPreintegratorBase::reset() {
    measurements_.clear();
    dT_ = 0;

    sigma_ = CovType::Zero(state_dim, state_dim);

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

} // namespace test_imu