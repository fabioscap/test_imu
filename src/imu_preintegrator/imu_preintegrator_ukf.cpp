#include "imu_preintegrator_ukf.h"

#include <srrg_geometry/geometry3d.h>

#define PHIidx 0 // phi
#define Vidx 3   // vel
#define Pidx 6   // pos
#define BAidx 9  // bias acc
#define BGidx 12 // bias gyro

#define NGidx 0  // noise gyro
#define NAidx 3  // noise acc
#define NBGidx 6 // bias random walk noise gyro
#define NBAidx 9 // bias random walk noise acc

// Those two values make the covariance singular
#define NBGIidx 12 // bias initial uncertainty gyro
#define NBAIidx 15 // bias initial uncertainty acc

using namespace test_imu;

void ImuPreintegratorUKF::preintegrate(const ImuMeasurement& m, float dt) {
  measurements_.push_back(m);

  // UKF: alias variables
  const core::Matrix3f& deltaR    = delta_incr_.get<0>().data();
  const core::Vector3f& deltaV    = delta_incr_.get<1>().data();
  const core::Vector3f& deltaP    = delta_incr_.get<2>().data();
  const core::Vector3f& bias_acc  = delta_incr_.get<3>().data();
  const core::Vector3f& bias_gyro = delta_incr_.get<4>().data();

  // correct the measurements
  core::Vector3f acc_c     = m.acceleration - bias_acc_;
  core::Vector3f ang_vel_c = m.angular_vel - bias_gyro_;

  // auxiliary variables
  core::Vector3f dtheta = ang_vel_c * dt;
  auto acc_skew         = core::geometry3d::skew(acc_c);
  auto dR               = core::geometry3d::expMapSO3(dtheta);
  core::fixRotation(dR);
  auto Jr = core::geometry3d::jacobianExpMapSO3(dtheta);

  /* bias correction jacobians */
  // Position
  dp_db_acc_ += dv_db_acc_ * dt - 0.5 * deltaR * dt * dt;
  dp_db_gyro_ += dv_db_gyro_ * dt - 0.5 * deltaR * acc_skew * dR_db_gyro_ * dt * dt;

  // Velocity
  dv_db_acc_ -= deltaR * dt;
  dv_db_gyro_ -= deltaR * acc_skew * dR_db_gyro_ * dt;

  // Rotation
  // remember that dR_j_j = I
  dR_db_gyro_ = dR.transpose() * dR_db_gyro_ - Jr * dt;

  scaling_.setIdentity(); // necessary?
  // it can also be precomputed if dt is constant
  scaling_.block<3, 3>(NGidx, NGidx) = Eigen::Matrix3f::Identity() / dt;
  scaling_.block<3, 3>(NAidx, NAidx) = Eigen::Matrix3f::Identity() / dt;
  // scaling_.block<3, 3>(NBGIidx, NBGIidx) = Eigen::Matrix3f::Identity() / dt;
  // scaling_.block<3, 3>(NBAIidx, NBAIidx) = Eigen::Matrix3f::Identity() / dt;
  scaling_.block<3, 3>(NBGidx, NBGidx) = Eigen::Matrix3f::Identity() * dt;
  scaling_.block<3, 3>(NBAidx, NBAidx) = Eigen::Matrix3f::Identity() * dt;

  /* estimate updates */
  sigma_joint_.setZero();
  sigma_joint_.block<state_dim, state_dim>(0, 0)                 = sigma_;
  sigma_joint_.block<noise_dim, noise_dim>(state_dim, state_dim) = scaling_ * sigma_noise_;

  // UKF: set input mean
  delta_incr_.get<5>().setData(m.acceleration);
  delta_incr_.get<6>().setData(m.angular_vel);

  UnscentedTransform::toUnscented(delta_incr_, sigma_joint_, spoints);

  // for each sigma points we do forward dynamics
  for (size_t i = 0; i < spoints.points.size(); ++i) {
    DeltaManifold& point                  = spoints.points.at(i);
    const core::Matrix3f& deltaR          = point.get<0>().data();
    const core::Vector3f& deltaV          = point.get<1>().data();
    const core::Vector3f& deltaP          = point.get<2>().data();
    const core::Vector3f& bias_acc        = point.get<3>().data();
    const core::Vector3f& bias_gyro       = point.get<4>().data();
    const core::Vector3f& acc_noisy_meas  = point.get<5>().data();
    const core::Vector3f& gyro_noisy_meas = point.get<6>().data();
    const core::Vector3f& noise_bacc      = point.get<7>().data();
    const core::Vector3f& noise_bgyro     = point.get<8>().data();

    // nominal values or UKF-generated???????
    auto acc_c     = acc_noisy_meas - bias_acc;
    auto ang_vel_c = gyro_noisy_meas - bias_gyro;

    // auxiliary variables
    core::Vector3f dtheta = ang_vel_c * dt;
    auto dR               = core::geometry3d::expMapSO3(dtheta);

    auto dva = deltaR * (acc_c) *dt;
    point.get<2>().setData(deltaP + (deltaV + 0.5 * dva) * dt);

    point.get<1>().setData(deltaV + dva);

    core::Matrix3f temp = deltaR * dR;
    core::fixRotation(temp);
    point.get<0>().setData(temp.eval());

    // bias evolution
    point.get<3>().setData(bias_acc + noise_bacc);
    point.get<4>().setData(bias_gyro + noise_bgyro);
  }

  // go back to normal parametrization
  UnscentedTransform::toMeanCov(spoints, delta_incr_, sigma_joint_);

  sigma_ = sigma_joint_.block<state_dim, state_dim>(0, 0);

  dT_ += dt;
}

void ImuPreintegratorUKF::reset() {
  ImuPreintegratorBase::reset();

  delta_incr_ = DeltaManifold();
}

void ImuPreintegratorUKF::getPrediction(const core::Isometry3f& Ti,
                                        const core::Vector3f& vi,
                                        core::Isometry3f& Tf,
                                        core::Vector3f& vf) {
  Tf.setIdentity();

  float T = measurements_.back().timestamp - measurements_.at(0).timestamp;

  vf = Ti.linear() * delta_v() + vi;

  Tf.linear()      = Ti.linear() * delta_R();
  Tf.translation() = Ti.linear() * delta_p() + Ti.translation() + T * vi;
}