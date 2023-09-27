#include "imu_preintegrator.h"

#include <srrg_geometry/geometry3d.h>

#define PHIidx 0 // phi
#define Vidx 3   // vel
#define Pidx 6   // pos
#define BAidx 9  // bias acc
#define BGidx 12 // bias gyro

#define NGidx 0    // noise gyro
#define NAidx 3    // noise acc
#define NBGidx 6   // bias random walk noise gyro
#define NBAidx 9   // bias random walk noise acc
#define NBGIidx 12 // bias initial uncertainty gyro
#define NBAIidx 15 // bias initial uncertainty acc

using namespace test_imu;

void ImuPreintegrator::preintegrate(const ImuMeasurement& m, float dt) {
  measurements_.push_back(m);

  // correct the measurements
  core::Vector3f acc_c     = m.acceleration - bias_acc_;
  core::Vector3f ang_vel_c = m.angular_vel - bias_gyro_;

  // auxiliary variables
  core::Vector3f dtheta = ang_vel_c * dt;
  auto acc_skew         = core::geometry3d::skew(acc_c);
  auto dR               = core::geometry3d::expMapSO3(dtheta);
  core::fixRotation(dR);
  auto Jr = core::geometry3d::jacobianExpMapSO3(dtheta);

  /* covariance matrix update through iterative first order propagation */
  // we do it before updating the estimates because it uses old delta_R

  // the biases jacobians are the same as the noise jacobians because
  // acc = acc_tilde - bias - noise
  // they are interchangeable

  // since preintegration needs biases, which are an estimation variable with associated covariance,
  // we also include the initial covariance of the estimate in the first order propagation

  // A_ ...
  A_.setIdentity(); // necessary?
  // dphi_dphi
  A_.block<3, 3>(PHIidx, PHIidx) = dR.transpose();
  // dphi_dv: 0, dphi_dp: 0
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
  A_.block<3, 3>(Pidx, Vidx) *= dt;
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
  B_.block<3, 3>(PHIidx, NBGIidx) = dt * Jr;
  // dphi_dbacc_noise_init: 0, dv_dgyro_noise: 0
  // dv_dacc_noise
  B_.block<3, 3>(Vidx, NAidx) = delta_R_ * dt;
  // dv_dbgyro_noise: 0, dv_dbacc_noise: 0, dv_dbgyro_noise_init: 0,
  // dv_dbacc_noise_init
  B_.block<3, 3>(Vidx, NBAIidx) = delta_R_ * dt;
  // dp_dgyro_noise: 0
  // dp_dacc_noise
  B_.block<3, 3>(Pidx, NAidx) = 0.5 * delta_R_ * dt * dt;
  // dp_dbgyro_noise: 0, dp_dbacc_noise: 0, dp_dbgyro_noise_init: 0,
  // dp_dbacc_noise_init
  B_.block<3, 3>(Pidx, NBAIidx) = 0.5 * delta_R_ * dt * dt;
  // dbgyro_dgyro_noise: 0, dbgyro_dacc_noise: 0
  // dbgyro_dbgyro_noise:
  B_.block<3, 3>(BGidx, NBGidx) = core::Matrix3f::Identity();
  // dbgyro_dbacc_noise: 0, dbgyro_dbgyro_noise_init: 0, dbgyro_dbacc_noise_init: 0,
  // dbacc_dgyro_noise: 0, dbacc_dacc_noise: 0, dbacc_dbgyro_noise: 0
  // dbacc_dbacc_noise
  B_.block<3, 3>(BAidx, NBAidx) = core::Matrix3f::Identity();
  // dbacc_dbgyro_noise_init: 0, dbacc_dbacc_noise_init: 0

  // TODO B is sparse and sigma_noise_ is block diagonal
  // it is possible to do this much quicker
  scaling_.setZero(); // necessary?
  // it can also be precomputed if dt is constant
  scaling_.block<3, 3>(NGidx, NGidx)     = Eigen::Matrix3f::Ones() / dt;
  scaling_.block<3, 3>(NAidx, NAidx)     = Eigen::Matrix3f::Ones() / dt;
  scaling_.block<3, 3>(NBGIidx, NBGIidx) = Eigen::Matrix3f::Ones() / dt;
  scaling_.block<3, 3>(NBAIidx, NBAIidx) = Eigen::Matrix3f::Ones() / dt;
  scaling_.block<3, 3>(NBGidx, NBGidx)   = Eigen::Matrix3f::Ones() * dt;
  scaling_.block<3, 3>(NBAidx, NBAidx)   = Eigen::Matrix3f::Ones() * dt;

  sigma_ = A_ * sigma_ * A_.transpose() + B_ * scaling_ * sigma_noise_ * B_.transpose();

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
  measurements_.clear();

  // initialize biases

  // set imu covariance matrix
}

void ImuPreintegrator::getPrediction(const core::Isometry3f& Ti,
                                     const core::Vector3f& vi,
                                     core::Isometry3f& Tf,
                                     core::Vector3f& vf) {
  Tf.setIdentity();

  float T = measurements_.back().timestamp - measurements_.at(0).timestamp;

  vf = Ti.linear() * delta_v_ + vi;

  Tf.linear()      = Ti.linear() * delta_R_;
  Tf.translation() = Ti.linear() * delta_p_ + Ti.translation() + T * vi;
}