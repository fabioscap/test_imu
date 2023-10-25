#pragma once

#include "common/common.h"

namespace test_imu {

  class ImuPreintegratorBase {
  public:
    // the data type used mainly for covariance propagation
    using Scalar = test_imu::Scalar;

    using Vector3 = core::Vector3_<Scalar>;
    using Matrix3 = core::Matrix3_<Scalar>;

    static constexpr int state_dim = 15;
    static constexpr int noise_dim = 12;

    static constexpr int PHIidx = 0;  // phi
    static constexpr int Vidx   = 3;  // vel
    static constexpr int Pidx   = 6;  // pos
    static constexpr int BAidx  = 9;  // bias acc
    static constexpr int BGidx  = 12; // bias gyro

    static constexpr int NGidx  = 0; // noise gyro
    static constexpr int NAidx  = 3; // noise acc
    static constexpr int NBGidx = 6; // bias random walk noise gyro
    static constexpr int NBAidx = 9; // bias random walk noise acc

    using CovType      = Eigen::Matrix<Scalar, -1, -1>;
    using CovNoiseType = Eigen::Matrix<Scalar, -1, -1>;

    virtual void preintegrate(const ImuMeasurement&, Scalar dt) = 0;
    virtual void reset();

    virtual const core::Matrix3f delta_R() const = 0;
    virtual const core::Vector3f delta_v() const = 0;
    virtual const core::Vector3f delta_p() const = 0;
    virtual const CovType& sigma() const         = 0;

    const void setNoiseGyroscope(const core::Vector3f& v);
    const void setNoiseAccelerometer(const core::Vector3f& v);
    const void setNoiseBiasGyroscope(const core::Vector3f& v);
    const void setNoiseBiasAccelerometer(const core::Vector3f& v);

    inline const core::Vector3f bias_acc() const {
      return bias_acc_.cast<float>();
    }
    inline const core::Vector3f bias_gyro() const {
      return bias_gyro_.cast<float>();
    }
    inline const core::Matrix3f dR_db_gyro() const {
      return dR_db_gyro_.cast<float>();
    }
    inline const core::Matrix3f dv_db_acc() const {
      return dv_db_acc_.cast<float>();
    }
    inline const core::Matrix3f dv_db_gyro() const {
      return dv_db_gyro_.cast<float>();
    }
    inline const core::Matrix3f dp_db_acc() const {
      return dp_db_acc_.cast<float>();
    }
    inline const core::Matrix3f dp_db_gyro() const {
      return dp_db_gyro_.cast<float>();
    }
    inline const float dT() const {
      return dT_;
    }
    inline const void setBiasAcc(const core::Vector3f& v) {
      bias_acc_ = v.cast<Scalar>();
    }
    inline const void setBiasGyro(const core::Vector3f& v) {
      bias_gyro_ = v.cast<Scalar>();
    }

  protected:
    std::vector<ImuMeasurement> measurements_;

    CovType sigma_ = CovType::Zero(state_dim, state_dim);

    // block diagonal
    // this should remain the same during operation
    CovNoiseType sigma_noise_ = 1e-3 * CovNoiseType::Identity(noise_dim, noise_dim);
    // we need discretize the noise covariances
    // see Optimal state estimation 8.1
    // or https://github.com/borglab/gtsam/blob/develop/doc/ImuFactor.pdf
    // preallocate matrix for scaling
    CovNoiseType scaling_ = CovNoiseType::Identity(noise_dim, noise_dim);

    // nominal values for the bias
    Vector3 bias_acc_  = Vector3::Zero();
    Vector3 bias_gyro_ = Vector3::Zero();

    // bias correction
    Matrix3 dR_db_gyro_ = Matrix3::Zero();
    Matrix3 dv_db_acc_  = Matrix3::Zero();
    Matrix3 dv_db_gyro_ = Matrix3::Zero();
    Matrix3 dp_db_acc_  = Matrix3::Zero();
    Matrix3 dp_db_gyro_ = Matrix3::Zero();

    Scalar dT_; // total preintegration time
  };
} // namespace test_imu