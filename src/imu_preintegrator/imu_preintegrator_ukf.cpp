#include "imu_preintegrator_ukf.h"

#include <srrg_geometry/geometry3d.h>

#include "common/manifold_impl.cpp"
#include "common/unscented_impl.cpp"

#include <random>

namespace test_imu {

  const core::Matrix3f ImuPreintegratorUKF::delta_R() const {
    return delta_incr_.get<0>().data().cast<float>();
  }
  const core::Vector3f ImuPreintegratorUKF::delta_p() const {
    return delta_incr_.get<2>().data().cast<float>();
  }
  const core::Vector3f ImuPreintegratorUKF::delta_v() const {
    return delta_incr_.get<1>().data().cast<float>();
  }

  void ImuPreintegratorUKF::preintegrate(const ImuMeasurement& m, Scalar dt) {
    measurements_.push_back(m);

    // UKF: alias variables
    const Matrix3& deltaR    = delta_incr_.get<0>().data();
    const Vector3& deltaV    = delta_incr_.get<1>().data();
    const Vector3& deltaP    = delta_incr_.get<2>().data();
    const Vector3& bias_acc  = delta_incr_.get<3>().data();
    const Vector3& bias_gyro = delta_incr_.get<4>().data();

    // correct the measurements
    Vector3 acc_c     = m.acceleration.cast<Scalar>() - bias_acc_;
    Vector3 ang_vel_c = m.angular_vel.cast<Scalar>() - bias_gyro_;

    // auxiliary variables
    Vector3 dtheta = ang_vel_c * dt;
    auto acc_skew  = core::geometry3d::skew(acc_c);
    auto dR        = core::geometry3d::expMapSO3(dtheta);
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
    scaling_.block<3, 3>(NGidx, NGidx) = Matrix3::Identity() / dt;
    scaling_.block<3, 3>(NAidx, NAidx) = Matrix3::Identity() / dt;
    // scaling_.block<3, 3>(NBGIidx, NBGIidx) = Matrix3::Identity() / dt;
    // scaling_.block<3, 3>(NBAIidx, NBAIidx) = Matrix3::Identity() / dt;
    scaling_.block<3, 3>(NBGidx, NBGidx) = Matrix3::Identity() * dt;
    scaling_.block<3, 3>(NBAidx, NBAidx) = Matrix3::Identity() * dt;

    /* estimate updates */
    sigma_joint_.setZero();
    sigma_joint_.block<state_dim, state_dim>(0, 0) = sigma_;
    sigma_joint_.block<input_dim, input_dim>(state_dim, state_dim) =
      (scaling_ * sigma_noise_).block<input_dim, input_dim>(0, 0);

    // UKF: set input mean
    delta_incr_.get<5>().setData(m.acceleration.cast<Scalar>());
    delta_incr_.get<6>().setData(m.angular_vel.cast<Scalar>());

    UnscentedTransform ut;
    ut.toUnscented(delta_incr_, sigma_joint_, spoints);

    // for each sigma points we do forward dynamics
    for (size_t i = 0; i < spoints.points.size(); ++i) {
      DeltaManifold& point           = spoints.points.at(i);
      const Matrix3& deltaR          = point.get<0>().data();
      const Vector3& deltaV          = point.get<1>().data();
      const Vector3& deltaP          = point.get<2>().data();
      const Vector3& bias_acc        = point.get<3>().data();
      const Vector3& bias_gyro       = point.get<4>().data();
      const Vector3& acc_noisy_meas  = point.get<5>().data();
      const Vector3& gyro_noisy_meas = point.get<6>().data();
      // const Vector3& noise_bacc      = point.get<7>().data();
      // const Vector3& noise_bgyro     = point.get<8>().data();

      // nominal values or UKF-generated???????
      auto acc_c     = acc_noisy_meas - bias_acc;
      auto ang_vel_c = gyro_noisy_meas - bias_gyro;

      // auxiliary variables
      Vector3 dtheta = ang_vel_c * dt;
      auto dR        = core::geometry3d::expMapSO3(dtheta);

      auto dva = deltaR * (acc_c) *dt;
      point.get<2>().setData(deltaP + (deltaV + 0.5 * dva) * dt);

      point.get<1>().setData(deltaV + dva);

      Matrix3 temp = deltaR * dR;
      core::fixRotation(temp);
      point.get<0>().setData(temp.eval());
    }

    // go back to normal parametrization
    ut.toMeanCov(spoints, delta_incr_, sigma_joint_);

    sigma_ = sigma_joint_.block<state_dim, state_dim>(0, 0);

    // add process noise
    sigma_.block<6, 6>(9, 9) += (scaling_ * sigma_noise_).block<6, 6>(6, 6);

    dT_ += dt;
  }

  void ImuPreintegratorUKF::reset() {
    ImuPreintegratorBase::reset();

    sigma_ = 1e-10 * CovType::Identity(state_dim, state_dim);

    delta_incr_ = DeltaManifold();
    sigma_joint_.setZero();
  }

  const core::Matrix3f ImuPreintegratorUKFSlim::delta_R() const {
    return delta_incr_.get<0>().data().cast<float>();
  }
  const core::Vector3f ImuPreintegratorUKFSlim::delta_p() const {
    return delta_incr_.get<2>().data().cast<float>();
  }
  const core::Vector3f ImuPreintegratorUKFSlim::delta_v() const {
    return delta_incr_.get<1>().data().cast<float>();
  }

  void ImuPreintegratorUKFSlim::preintegrate(const ImuMeasurement& m, Scalar dt) {
    measurements_.push_back(m);

    scaling_.setIdentity(); // necessary?
    // it can also be precomputed if dt is constant
    scaling_.block<3, 3>(NGidx, NGidx) = Matrix3::Identity() / dt;
    scaling_.block<3, 3>(NAidx, NAidx) = Matrix3::Identity() / dt;

    /* estimate updates */
    // sigma_joint_.setZero();
    // sigma_joint_.block<state_dim, state_dim>(0, 0) = sigma_;
    sigma_joint_.block<input_dim, input_dim>(state_dim, state_dim) =
      (scaling_ * sigma_noise_).block<input_dim, input_dim>(0, 0).array().sqrt();
    sigma_joint_.block<input_dim, state_dim>(state_dim, 0).setZero();
    sigma_joint_.block<state_dim, input_dim>(0, state_dim).setZero();
    /*     std::cout << "t: " << m.timestamp << "\n";
        std::cout << "sigma before propagation:\n" << sigma_joint_ << "\n";
     */
    // UKF: set input mean
    delta_incr_.get<3>().setData(m.acceleration.cast<Scalar>());
    delta_incr_.get<4>().setData(m.angular_vel.cast<Scalar>());

    UnscentedTransform ut;
    ut.toSqrtUnscented(delta_incr_, sigma_joint_, spoints);

    // for each sigma points we do forward dynamics
    for (size_t i = 0; i < spoints.points.size(); ++i) {
      DeltaManifold& point           = spoints.points.at(i);
      const Matrix3& deltaR          = point.get<0>().data();
      const Vector3& deltaV          = point.get<1>().data();
      const Vector3& deltaP          = point.get<2>().data();
      const Vector3& acc_noisy_meas  = point.get<3>().data();
      const Vector3& gyro_noisy_meas = point.get<4>().data();

      auto acc_c     = acc_noisy_meas - bias_acc_;
      auto ang_vel_c = gyro_noisy_meas - bias_gyro_;

      // auxiliary variables
      Vector3 dtheta = ang_vel_c * dt;
      auto dR        = core::geometry3d::expMapSO3(dtheta);

      auto dva = deltaR * (acc_c) *dt;
      point.get<2>().setData(deltaP + (deltaV + 0.5 * dva) * dt);

      point.get<1>().setData(deltaV + dva);

      Matrix3 temp = deltaR * dR;
      core::fixRotation(temp);
      point.get<0>().setData(temp.eval());
    }

    ut.toMeanSqrtCov(spoints, delta_incr_, sigma_joint_);

    sigma_ = (sigma_joint_ * sigma_joint_.transpose()).block<state_dim, state_dim>(0, 0);

    dT_ += dt;
  }

  void ImuPreintegratorUKFSlim::reset() {
    ImuPreintegratorBase::reset();

    sigma_ = 1e-10 * CovType::Identity(state_dim, state_dim);

    sigma_joint_.setIdentity();
    sigma_joint_ *= 1e-5;
    delta_incr_ = DeltaManifold();
  }

  ImuPreintegratorMC::Vector3 ImuPreintegratorMC::sample_noise(const Vector3& cov_diag) {
    // Create a random number generator
    std::default_random_engine generator;

    // Create normal distributions for each dimension
    std::normal_distribution<Scalar> dist_x(0.0, std::sqrt(cov_diag(0)));
    std::normal_distribution<Scalar> dist_y(0.0, std::sqrt(cov_diag(1)));
    std::normal_distribution<Scalar> dist_z(0.0, std::sqrt(cov_diag(2)));

    Vector3 noise;
    noise << dist_x(generator), dist_y(generator), dist_z(generator);

    return noise;
  }

  void ImuPreintegratorMC::compute_sigma_from_particles() {
    sigma_ = core::MatrixN_<Scalar, 9>::Zero();
    for (size_t i = 0; i < particles_.size(); ++i) {
      auto err = particles_.at(i).boxminus(delta_incr_);
      sigma_ += err * err.transpose();
    }
    sigma_ /= particles_.size();
  }

  void ImuPreintegratorMC::preintegrate(const ImuMeasurement& m, Scalar dt) {
    measurements_.push_back(m);

    // UKF: alias variables
    const Matrix3& deltaR = delta_incr_.get<0>().data();
    const Vector3& deltaV = delta_incr_.get<1>().data();
    const Vector3& deltaP = delta_incr_.get<2>().data();

    // correct the measurements
    Vector3 acc_c     = m.acceleration.cast<Scalar>() - bias_acc_;
    Vector3 ang_vel_c = m.angular_vel.cast<Scalar>() - bias_gyro_;

    // auxiliary variables
    Vector3 dtheta = ang_vel_c * dt;
    auto dR        = core::geometry3d::expMapSO3(dtheta);

    auto dva = deltaR * (acc_c) *dt;
    delta_incr_.get<2>().setData(deltaP + (deltaV + 0.5 * dva) * dt);

    delta_incr_.get<1>().setData(deltaV + dva);

    Matrix3 temp = deltaR * dR;
    core::fixRotation(temp);
    delta_incr_.get<0>().setData(temp.eval());

    for (size_t i = 0; i < particles_.size(); ++i) {
      const Matrix3& deltaR = particles_.at(i).get<0>().data();
      const Vector3& deltaV = particles_.at(i).get<1>().data();
      const Vector3& deltaP = particles_.at(i).get<2>().data();
      // perturb inputs
      const Vector3 noisy_acc     = acc_c + sample_noise(sigma_noise_.diagonal().head(3) / dt);
      const Vector3 noisy_ang_vel = ang_vel_c + sample_noise(sigma_noise_.diagonal().tail(3) / dt);

      // auxiliary variables
      Vector3 dtheta = noisy_ang_vel * dt;
      auto dR        = core::geometry3d::expMapSO3(dtheta);

      auto dva = deltaR * (noisy_acc) *dt;
      particles_.at(i).get<2>().setData(deltaP + (deltaV + 0.5 * dva) * dt);

      particles_.at(i).get<1>().setData(deltaV + dva);

      Matrix3 temp = deltaR * dR;
      core::fixRotation(temp);
      particles_.at(i).get<0>().setData(temp.eval());
    }

    dT_ += dt;

    compute_sigma_from_particles();
  }

  void ImuPreintegratorMC::reset() {
    ImuPreintegratorBase::reset();
    sigma_.setZero();
    delta_incr_ = DeltaManifold();
  }

  const core::Matrix3f ImuPreintegratorMC::delta_R() const {
    return delta_incr_.get<0>().data().cast<float>();
  }
  const core::Vector3f ImuPreintegratorMC::delta_p() const {
    return delta_incr_.get<2>().data().cast<float>();
  }
  const core::Vector3f ImuPreintegratorMC::delta_v() const {
    return delta_incr_.get<1>().data().cast<float>();
  }

} // namespace test_imu