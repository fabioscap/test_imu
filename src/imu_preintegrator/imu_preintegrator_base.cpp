#include "imu_preintegrator_base.h"

namespace test_imu {
  void ImuPreintegratorBase::preintegrate(const ImuMeasurement& m, Scalar dt) {
    measurements_.push_back(m);
    dT_ += dt;
    preintegrate_(m, dt);
  }

  void ImuPreintegratorBase::reset() {
    measurements_.clear();
    dT_ = 0;
    reset_();
  }

  void ImuPreintegratorBase::f(Matrix3& delta_R,
                               Vector3& delta_v,
                               Vector3& delta_p,
                               const Matrix3& dR,
                               const Vector3& acc,
                               float dt) {
    Vector3 dva = delta_R * (acc) *dt;
    delta_p += (delta_v + 0.5 * dva) * dt;

    // update the velocity
    delta_v += dva;

    // update orientation
    delta_R *= dR;
    core::fixRotation(delta_R);
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

  const void ImuPreintegratorBase::setNoiseGyroscope(const core::Vector3f& v) {
    gyro_noise_ = (v.cast<Scalar>());
  }
  const void ImuPreintegratorBase::setNoiseAccelerometer(const core::Vector3f& v) {
    acc_noise_ = (v.cast<Scalar>());
  }
  const void ImuPreintegratorBase::setNoiseBiasGyroscope(const core::Vector3f& v) {
    bias_gyro_noise_ = (v.cast<Scalar>());
  }
  const void ImuPreintegratorBase::setNoiseBiasAccelerometer(const core::Vector3f& v) {
    bias_acc_noise_ = (v.cast<Scalar>());
  }

  const void ImuPreintegratorBase::setBiasAcc(const core::Vector3f& v) {
    bias_acc_ = v.cast<Scalar>();
  }
  const void ImuPreintegratorBase::setBiasGyro(const core::Vector3f& v) {
    bias_gyro_ = v.cast<Scalar>();
  }

  const float ImuPreintegratorBase::dT() const {
    return dT_;
  }

  const core::Vector3f ImuPreintegratorBase::bias_gyro() const {
    return bias_gyro_.cast<float>();
  }

  const core::Vector3f ImuPreintegratorBase::bias_acc() const {
    return bias_acc_.cast<float>();
  }

  std::vector<ImuMeasurement> ImuPreintegratorBase::measurements() {
    return measurements_;
  }

  void BiasJacobians::update(const Matrix3& deltaR,
                             const Matrix3& acc_skew,
                             const Matrix3& dR,
                             const Matrix3& Jr,
                             Scalar dt) {
    dp_db_acc += dv_db_acc * dt - 0.5 * deltaR * dt * dt;
    dp_db_gyro += dv_db_gyro * dt - 0.5 * deltaR * acc_skew * dR_db_gyro * dt * dt;

    // Velocity
    dv_db_acc -= deltaR * dt;
    dv_db_gyro -= deltaR * acc_skew * dR_db_gyro * dt;

    // Rotation
    // remember that dR_j_j = I
    dR_db_gyro = dR.transpose() * dR_db_gyro - Jr * dt;
  }
} // namespace test_imu