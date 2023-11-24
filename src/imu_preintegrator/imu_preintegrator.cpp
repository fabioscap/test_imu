#include "imu_preintegrator.h"

#include <srrg_geometry/geometry3d.h>

#define PHIidx 0
#define Vidx 3
#define Pidx 6
#define BAidx 9
#define BGidx 12

#define NGidx 0
#define NAidx 3
#define NBGidx 6
#define NBAidx 9

namespace test_imu {
  const BiasJacobians* ImuPreintegrator::biasJacobians() const {
    return &bias_J_;
  }

  void ImuPreintegrator::preintegrate_(const ImuMeasurement& m, Scalar dt) {
    // correct the measurements
    Vector3 acc_c     = m.acceleration.cast<Scalar>() - bias_acc_;
    Vector3 ang_vel_c = m.angular_vel.cast<Scalar>() - bias_gyro_;

    // auxiliary variables
    Vector3 dtheta = ang_vel_c * dt;
    auto acc_skew  = core::geometry3d::skew(acc_c);
    auto dR        = core::geometry3d::expMapSO3(dtheta);
    core::fixRotation(dR);
    auto Jr = core::geometry3d::jacobianExpMapSO3(dtheta);

    /* covariance matrix update through iterative first order propagation */
    // we do it before updating the estimates because it uses old delta_R

    // the biases jacobians are the same as the noise jacobians because
    // acc = acc_tilde - bias - noise
    // they are interchangeable

    // A_ ...
    // A_.setIdentity();
    // dphi_dphi
    A_.block<3, 3>(PHIidx, PHIidx) = dR.transpose();
    // dphi_dv: 0, dphi_dp: 0, dphi_dbacc: 0
    // dphi_dbgyro
    A_.block<3, 3>(PHIidx, BGidx) = dt * Jr;
    // dv_dphi
    A_.block<3, 3>(Vidx, PHIidx) = -delta_R_ * acc_skew * dt;
    // dv_dv: I, dv_dp: 0, dv_dbgyro: 0,
    // dv_dbacc
    A_.block<3, 3>(Vidx, BAidx) = delta_R_ * dt;
    // dp_dphi
    A_.block<3, 3>(Pidx, PHIidx) = -0.5 * delta_R_ * acc_skew * dt * dt;
    // dp_dv
    A_.block<3, 3>(Pidx, Vidx) = Matrix3::Identity() * dt;
    // dp_dp: I , dp_dbgyro: 0
    // dp_dbacc
    A_.block<3, 3>(Pidx, BAidx) = 0.5 * delta_R_ * dt * dt;
    // dbgyro_dphi: 0, dbgyro_dv: 0, dbgyro_dp: 0, dbgyro_dbgyro: I, dbgyro_dbacc: 0, dbacc_dphi: 0,
    // dbacc_dv: 0, dbacc_dp: 0, dbacc_dbgyro: 0, dbacc_dbacc: I

    // B_ ...
    // B_.setZero();
    // dphi_dgyro_noise
    B_.block<3, 3>(PHIidx, NGidx) = dt * Jr;
    // dphi_dacc_noise: 0,dphi_dbgyro_noise: 0, dphi_dbacc_noise: 0
    // dv_dacc_noise:
    B_.block<3, 3>(Vidx, NAidx) = delta_R_ * dt;
    // dv_dbgyro_noise: 0, dv_dbacc_noise: 0, dp_dgyro_noise: 0
    // dp_dacc_noise:
    B_.block<3, 3>(Pidx, NAidx) = 0.5 * delta_R_ * dt * dt;
    // dp_dbgyro_noise: 0, dp_dbacc_noise: 0,
    // dbgyro_dgyro_noise: 0, dbgyro_dacc_noise: 0,
    // dbgyro_dbgyro_noise:
    B_.block<3, 3>(BGidx, NBGidx) = Matrix3::Identity();
    // dbgyro_dbacc_noise: 0, dbacc_dgyro_noise: 0, dbacc_dacc_noise: 0, dbacc_dbgyro_noise: 0,
    // dbacc_dbacc_noise:
    B_.block<3, 3>(BAidx, NBAidx) = Matrix3::Identity();

    sigma_noise_diag_ << gyro_noise_ / dt, acc_noise_ / dt, bias_gyro_noise_ * dt,
      bias_acc_noise_ * dt;

    sigma_ = A_ * sigma_ * A_.transpose() + B_ * sigma_noise_diag_.asDiagonal() * B_.transpose();

    bias_J_.update(delta_R_, acc_skew, dR, Jr, dt);
    f(delta_R_, delta_v_, delta_p_, dR, acc_c, dt);
  }

  void ImuPreintegrator::reset_() {
    sigma_ = 1e-10 * CovType::Identity();

    delta_p_ = Vector3::Zero();
    delta_R_ = Matrix3::Identity();
    delta_v_ = Vector3::Zero();

    // allocate matrices for noise propagation
    A_.setIdentity();
    B_.setZero();

    bias_J_ = BiasJacobians();
  }

  void ImuPreintegratorSlim::preintegrate_(const ImuMeasurement& m, Scalar dt) {
    // correct the measurements
    Vector3 acc_c     = m.acceleration.cast<Scalar>() - bias_acc_;
    Vector3 ang_vel_c = m.angular_vel.cast<Scalar>() - bias_gyro_;

    // auxiliary variables
    Vector3 dtheta = ang_vel_c * dt;
    auto acc_skew  = core::geometry3d::skew(acc_c);
    auto dR        = core::geometry3d::expMapSO3(dtheta);
    core::fixRotation(dR);
    auto Jr = core::geometry3d::jacobianExpMapSO3(dtheta);

    // A_ ...
    A_.setIdentity(); // necessary?
    A_.block<3, 3>(PHIidx, PHIidx) = dR.transpose();
    A_.block<3, 3>(Vidx, PHIidx)   = -delta_R_ * acc_skew * dt;
    A_.block<3, 3>(Pidx, PHIidx)   = -0.5 * delta_R_ * acc_skew * dt * dt;
    A_.block<3, 3>(Pidx, Vidx)     = Matrix3::Identity() * dt;

    // B_ ...
    B_.setZero(); // necessary?
    B_.block<3, 3>(PHIidx, NGidx) = dt * Jr;
    B_.block<3, 3>(Vidx, NAidx)   = delta_R_ * dt;
    B_.block<3, 3>(Pidx, NAidx)   = 0.5 * delta_R_ * dt * dt;

    sigma_noise_diag_ << gyro_noise_ / dt, acc_noise_ / dt;

    sigma_ = A_ * sigma_ * A_.transpose() + B_ * sigma_noise_diag_.asDiagonal() * B_.transpose();

    f(delta_R_, delta_v_, delta_p_, dR, acc_c, dt);
  }

  void ImuPreintegratorSlim::reset_() {
    sigma_ = 1e-10 * CovType::Identity(state_dim, state_dim);

    delta_p_ = Vector3::Zero();
    delta_R_ = Matrix3::Identity();
    delta_v_ = Vector3::Zero();

    // allocate matrices for noise propagation
    A_.setIdentity();
    B_.setZero();
  }

} // namespace test_imu