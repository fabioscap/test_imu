#include "imu_preintegrator_slim.h"

#include <srrg_geometry/geometry3d.h>

#define PHIidx 0 // phi
#define Vidx 3   // vel
#define Pidx 6   // pos

#define NGidx 0  // noise gyro
#define NAidx 3  // noise acc
#define NBGidx 6 // bias random walk noise gyro
#define NBAidx 9 // bias random walk noise acc

namespace test_imu {
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

    /* covariance matrix update through iterative first order propagation */
    // we do it before updating the estimates because it uses old delta_R

    // the biases jacobians are the same as the noise jacobians because
    // acc = acc_tilde - bias - noise
    // they are interchangeable

    // since preintegration needs biases, which are an estimation variable with associated
    // covariance, we also include the initial covariance of the estimate in the first order
    // propagation

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
    /*     scaling_.block<3, 3>(NBGidx, NBGidx) = Matrix3::Identity() * dt;
        scaling_.block<3, 3>(NBAidx, NBAidx) = Matrix3::Identity() * dt;
     */

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

    delta_p_ = Vector3::Zero();
    delta_R_ = Matrix3::Identity();
    delta_v_ = Vector3::Zero();

    // allocate matrices for noise propagation
    A_.setIdentity();
    B_.setIdentity();
  }

  void ImuPreintegratorSlim::getPrediction(const core::Isometry3f& Ti,
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