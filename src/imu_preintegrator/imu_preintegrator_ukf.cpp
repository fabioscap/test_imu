#include "imu_preintegrator_ukf.h"

#include <srrg_geometry/geometry3d.h>

#include "common/manifold_impl.cpp"
#include "common/unscented_impl.cpp"

#include <random>

namespace test_imu {

  const core::Matrix3f ImuPreintegratorUKF::delta_R() {
    return delta_incr_.get<0>().data().cast<float>();
  }
  const core::Vector3f ImuPreintegratorUKF::delta_p() {
    return delta_incr_.get<2>().data().cast<float>();
  }
  const core::Vector3f ImuPreintegratorUKF::delta_v() {
    return delta_incr_.get<1>().data().cast<float>();
  }
  const core::MatrixX_<float> ImuPreintegratorUKF::sigma() {
    return sigma_joint_.block<state_dim, state_dim>(0, 0).cast<float>();
  }

  void ImuPreintegratorUKF::preintegrate_(const ImuMeasurement& m, Scalar dt) {
    // UKF: alias variables
    const Matrix3& deltaR = delta_incr_.get<0>().data();

    // correct the measurements
    Vector3 acc_c     = m.acceleration.cast<Scalar>() - bias_acc_;
    Vector3 ang_vel_c = m.angular_vel.cast<Scalar>() - bias_gyro_;

    // auxiliary variables
    Vector3 dtheta   = ang_vel_c * dt;
    Matrix3 acc_skew = core::geometry3d::skew(acc_c);
    Matrix3 dR       = core::geometry3d::expMapSO3(dtheta);
    core::fixRotation(dR);
    Matrix3 Jr = core::geometry3d::jacobianExpMapSO3(dtheta);

    /* bias correction jacobians */
    bias_J_.update(deltaR, acc_skew, dR, Jr, dt);

    /* generation of sigma points */

    input_noise_diag_ << acc_noise_ / dt, gyro_noise_ / dt;
    bias_noise_diag_ << bias_acc_noise_ * dt, bias_gyro_noise_ * dt;

    delta_incr_.get<5>().data() = (m.acceleration.cast<Scalar>());
    delta_incr_.get<6>().data() = (m.angular_vel.cast<Scalar>());

    // zero out off diagonal blocks: x_t and u_t are independent
    sigma_joint_.block<input_dim, state_dim>(state_dim, 0).setZero();
    sigma_joint_.block<state_dim, input_dim>(0, state_dim).setZero();

    // set input noise
    sigma_joint_.block<input_dim, input_dim>(state_dim, state_dim) = input_noise_diag_.asDiagonal();

    ut_.toUnscented(delta_incr_, sigma_joint_, spoints);

    // for each sigma points we do forward dynamics
    for (size_t i = 0; i < spoints.points.size(); ++i) {
      DeltaManifold& point           = spoints.points.at(i);
      Matrix3& deltaR                = point.get<0>().data();
      Vector3& deltaV                = point.get<1>().data();
      Vector3& deltaP                = point.get<2>().data();
      Vector3& bias_acc              = point.get<3>().data();
      Vector3& bias_gyro             = point.get<4>().data();
      const Vector3& acc_noisy_meas  = point.get<5>().data();
      const Vector3& gyro_noisy_meas = point.get<6>().data();

      Vector3 acc_c     = acc_noisy_meas - bias_acc;
      Vector3 ang_vel_c = gyro_noisy_meas - bias_gyro;

      // auxiliary variables
      Vector3 dtheta = ang_vel_c * dt;
      Matrix3 dR     = core::geometry3d::expMapSO3(dtheta);

      // side-effect on sigma points
      f(deltaR, deltaV, deltaP, dR, acc_c, dt);
    }

    // go back to normal parametrization
    ut_.toMeanCov(spoints, delta_incr_, sigma_joint_);

    // add process noise
    sigma_joint_.block<6, 6>(9, 9) += bias_noise_diag_.asDiagonal();
  }

  void ImuPreintegratorUKF::reset_() {
    delta_incr_  = DeltaManifold();
    sigma_joint_ = 1e-10 * CovJType::Identity();

    spoints = SigmaPoints<DeltaManifold>();
  }

  const core::Matrix3f ImuPreintegratorUKFSlim::delta_R() {
    if (!is_updated_)
      toMeanCov_();
    return delta_incr_.get<0>().data().cast<float>();
  }
  const core::Vector3f ImuPreintegratorUKFSlim::delta_p() {
    if (!is_updated_)
      toMeanCov_();
    return delta_incr_.get<2>().data().cast<float>();
  }
  const core::Vector3f ImuPreintegratorUKFSlim::delta_v() {
    if (!is_updated_)
      toMeanCov_();
    return delta_incr_.get<1>().data().cast<float>();
  }
  const core::MatrixX_<float> ImuPreintegratorUKFSlim::sigma() {
    if (!is_updated_)
      toMeanCov_();
    return sigma_joint_.block<state_dim, state_dim>(0, 0).cast<float>();
  }

  void ImuPreintegratorUKFSlim::preintegrate_(const ImuMeasurement& m, Scalar dt) {
    // UKF: alias variables
    const Matrix3& deltaR = delta_incr_.get<0>().data();

    // correct the measurements
    Vector3 acc_c     = m.acceleration.cast<Scalar>() - bias_acc_;
    Vector3 ang_vel_c = m.angular_vel.cast<Scalar>() - bias_gyro_;

    // auxiliary variables
    Vector3 dtheta   = ang_vel_c * dt;
    Matrix3 acc_skew = core::geometry3d::skew(acc_c);
    Matrix3 dR       = core::geometry3d::expMapSO3(dtheta);
    core::fixRotation(dR);
    Matrix3 Jr = core::geometry3d::jacobianExpMapSO3(dtheta);

    /* bias correction jacobians */
    bias_J_.update(deltaR, acc_skew, dR, Jr, dt);

    if (measurements().size() == 1 || true)
      toUnscented_(dt);

    // for each sigma points we do forward dynamics
    for (size_t i = 0; i < spoints.points.size(); ++i) {
      DeltaManifold& point           = spoints.points.at(i);
      Matrix3& deltaR                = point.get<0>().data();
      Vector3& deltaV                = point.get<1>().data();
      Vector3& deltaP                = point.get<2>().data();
      const Vector3& acc_noisy_meas  = point.get<3>().data() + m.acceleration.cast<Scalar>();
      const Vector3& gyro_noisy_meas = point.get<4>().data() + m.angular_vel.cast<Scalar>();

      Vector3 acc_c     = acc_noisy_meas - bias_acc_;
      Vector3 ang_vel_c = gyro_noisy_meas - bias_gyro_;

      // auxiliary variables
      Vector3 dtheta = ang_vel_c * dt;
      Matrix3 dR     = core::geometry3d::expMapSO3(dtheta);

      // side-effect on sigma points
      f(deltaR, deltaV, deltaP, dR, acc_c, dt);
      // once we transform sigma points, mean and covariance are not updated
      is_updated_ = false;
    }
    toMeanCov_();
  }

  void ImuPreintegratorUKFSlim::reset_() {
    delta_incr_  = DeltaManifold();
    sigma_joint_ = 1e-10 * CovJType::Identity();

    spoints = SigmaPoints<DeltaManifold>();

    bias_J_     = BiasJacobians();
    is_updated_ = true;
  }

  void ImuPreintegratorUKFSlim::toUnscented_(float dt) {
    /* generation of sigma points */
    input_noise_diag_ << acc_noise_ / dt, gyro_noise_ / dt;

    // zero out off diagonal blocks: x_t and u_t are independent
    sigma_joint_.block<input_dim, state_dim>(state_dim, 0).setZero();
    sigma_joint_.block<state_dim, input_dim>(0, state_dim).setZero();

    // set input noise
    sigma_joint_.block<input_dim, input_dim>(state_dim, state_dim) = input_noise_diag_.asDiagonal();

    ut_.toUnscented(delta_incr_, sigma_joint_, spoints);
  }

  void ImuPreintegratorUKFSlim::toMeanCov_() {
    // go back to normal parametrization
    ut_.toMeanCov(spoints, delta_incr_, sigma_joint_);
    is_updated_ = true;
  }

} // namespace test_imu