#include "imu_preintegrator.h"

#include <srrg_geometry/geometry3d.h>

// Those two values make the covariance singular
// #define NBGIidx 12 // bias initial uncertainty gyro
// #define NBAIidx 15 // bias initial uncertainty acc

namespace test_imu {

  void ImuPreintegrator::preintegrate(const ImuMeasurement& m, Scalar dt) {
    measurements_.push_back(m);

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

    // since preintegration needs biases, which are an estimation variable with associated
    // covariance, we also include the initial covariance of the estimate in the first order
    // propagation // NO

    // A_ ...
    A_.setIdentity(); // necessary?
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
    B_.setZero(); // necessary?
    // dphi_dgyro_noise
    B_.block<3, 3>(PHIidx, NGidx) = dt * Jr;
    // dphi_dacc_noise: 0,dphi_dbgyro_noise: 0, dphi_dbacc_noise: 0
    // dphi_dbgyro_noise_init
    // B_.block<3, 3>(PHIidx, NBGIidx) = dt * Jr;
    // dphi_dbacc_noise_init: 0, dv_dgyro_noise: 0
    // dv_dacc_noise
    B_.block<3, 3>(Vidx, NAidx) = delta_R_ * dt;
    // dv_dbgyro_noise: 0, dv_dbacc_noise: 0, dv_dbgyro_noise_init: 0,
    // dv_dbacc_noise_init
    // B_.block<3, 3>(Vidx, NBAIidx) = delta_R_ * dt;
    // dp_dgyro_noise: 0
    // dp_dacc_noise
    B_.block<3, 3>(Pidx, NAidx) = 0.5 * delta_R_ * dt * dt;
    // dp_dbgyro_noise: 0, dp_dbacc_noise: 0, dp_dbgyro_noise_init: 0,
    // dp_dbacc_noise_init
    // B_.block<3, 3>(Pidx, NBAIidx) = 0.5 * delta_R_ * dt * dt;
    // dbgyro_dgyro_noise: 0, dbgyro_dacc_noise: 0
    // dbgyro_dbgyro_noise:
    B_.block<3, 3>(BGidx, NBGidx) = Matrix3::Identity();
    // dbgyro_dbacc_noise: 0, dbgyro_dbgyro_noise_init: 0, dbgyro_dbacc_noise_init: 0,
    // dbacc_dgyro_noise: 0, dbacc_dacc_noise: 0, dbacc_dbgyro_noise: 0
    // dbacc_dbacc_noise
    B_.block<3, 3>(BAidx, NBAidx) = Matrix3::Identity();
    // dbacc_dbgyro_noise_init: 0, dbacc_dbacc_noise_init: 0

    // TODO B is sparse and sigma_noise_ is block diagonal
    // it is possible to do this much quicker
    scaling_.setIdentity(); // necessary?
    // it can also be precomputed if dt is constant
    scaling_.block<3, 3>(NGidx, NGidx) = Matrix3::Identity() / dt;
    scaling_.block<3, 3>(NAidx, NAidx) = Matrix3::Identity() / dt;
    // scaling_.block<3, 3>(NBGIidx, NBGIidx) = Matrix3::Identity() / dt;
    // scaling_.block<3, 3>(NBAIidx, NBAIidx) = Matrix3::Identity() / dt;
    scaling_.block<3, 3>(NBGidx, NBGidx) = Matrix3::Identity() * dt;
    scaling_.block<3, 3>(NBAidx, NBAidx) = Matrix3::Identity() * dt;

    /*   std::cout << "B:\n" << B_ << std::endl;
      std::cout << "scaling_* sigma_noise_\n" << scaling_ * sigma_noise_ << "\n";
      std::cout << "Det\n"
                << (scaling_ * sigma_noise_).determinant() << "->"
                << (B_ * scaling_ * sigma_noise_ * B_.transpose()).determinant() << "\n";
     */
    sigma_ = A_ * sigma_ * A_.transpose() + B_ * scaling_ * sigma_noise_ * B_.transpose();
    // sigma_ = CovType::Identity(state_dim, state_dim);
    /* bias correction jacobians */
    // Position
    dp_db_acc_ += dv_db_acc_ * dt - 0.5 * delta_R_ * dt * dt;
    dp_db_gyro_ += dv_db_gyro_ * dt - 0.5 * delta_R_ * acc_skew * dR_db_gyro_ * dt * dt;

    // Velocity
    dv_db_acc_ -= delta_R_ * dt;
    dv_db_gyro_ -= delta_R_ * acc_skew * dR_db_gyro_ * dt;

    // Rotation
    // remember that dR_j_j = I
    dR_db_gyro_ = dR.transpose() * dR_db_gyro_ - Jr * dt;

    /* estimate updates */
    // update the position
    auto dva = delta_R_ * (acc_c) *dt;
    delta_p_ += (delta_v_ + 0.5 * dva) * dt;

    // update the velocity
    delta_v_ += dva;

    // update orientation
    delta_R_ *= dR;
    core::fixRotation(delta_R_);

    dT_ += dt;
  }

  void ImuPreintegrator::reset() {
    ImuPreintegratorBase::reset();
    sigma_ = 1e-10 * CovType::Identity(state_dim, state_dim);

    delta_p_ = Vector3::Zero();
    delta_R_ = Matrix3::Identity();
    delta_v_ = Vector3::Zero();

    // allocate matrices for noise propagation
    A_.setIdentity();
    B_.setZero();
  }

  void ImuPreintegratorSlim::preintegrate(const ImuMeasurement& m, Scalar dt) {
    measurements_.push_back(m);

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

    // it is possible to do this much quicker
    scaling_.setIdentity(); // necessary?
    // it can also be precomputed if dt is constant
    scaling_.block<3, 3>(NGidx, NGidx) = Matrix3::Identity() / dt;
    scaling_.block<3, 3>(NAidx, NAidx) = Matrix3::Identity() / dt;

    sigma_ = A_ * sigma_ * A_.transpose() + B_ * scaling_ * sigma_noise_ * B_.transpose();

    /* estimate updates */
    // update the position
    auto dva = delta_R_ * (acc_c) *dt;
    delta_p_ += (delta_v_ + 0.5 * dva) * dt;

    // update the velocity
    delta_v_ += dva;

    // update orientation
    delta_R_ *= dR;
    core::fixRotation(delta_R_);

    dT_ += dt;
  }
  void ImuPreintegratorSlim::reset() {
    ImuPreintegratorBase::reset();

    sigma_ =

      delta_p_ = Vector3::Zero();
    delta_R_   = Matrix3::Identity();
    delta_v_   = Vector3::Zero();

    // allocate matrices for noise propagation
    A_.setIdentity();
    B_.setZero();
  }

} // namespace test_imu