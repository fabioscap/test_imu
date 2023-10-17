#pragma once

#include <common/common.h>

namespace test_imu {

  class ImuPreintegratorBase {
  public:
    static constexpr int state_dim = 15;
    static constexpr int noise_dim = 12;

    using CovType      = Eigen::MatrixXf;
    using CovNoiseType = Eigen::MatrixXf;

    virtual void preintegrate(const ImuMeasurement&, float dt) = 0;
    virtual void reset();

    virtual const core::Matrix3f& delta_R() const = 0;
    virtual const core::Vector3f& delta_v() const = 0;
    virtual const core::Vector3f& delta_p() const = 0;
    virtual const CovType& sigma() const          = 0;

    const void setNoiseGyroscope(const core::Vector3f& v);
    const void setNoiseAccelerometer(const core::Vector3f& v);
    const void setNoiseBiasGyroscope(const core::Vector3f& v);
    const void setNoiseBiasAccelerometer(const core::Vector3f& v);

    inline const core::Vector3f& bias_acc() const {
      return bias_acc_;
    }
    inline const core::Vector3f& bias_gyro() const {
      return bias_gyro_;
    }
    inline const core::Matrix3f& dR_db_gyro() const {
      return dR_db_gyro_;
    }
    inline const core::Matrix3f& dv_db_acc() const {
      return dv_db_acc_;
    }
    inline const core::Matrix3f& dv_db_gyro() const {
      return dv_db_gyro_;
    }
    inline const core::Matrix3f& dp_db_acc() const {
      return dp_db_acc_;
    }
    inline const core::Matrix3f& dp_db_gyro() const {
      return dp_db_gyro_;
    }
    inline const float dT() const {
      return dT_;
    }
    inline const void setBiasAcc(const core::Vector3f& v) {
      bias_acc_ = v;
    }
    inline const void setBiasGyro(const core::Vector3f& v) {
      bias_gyro_ = v;
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
    core::Vector3f bias_acc_  = core::Vector3f::Zero();
    core::Vector3f bias_gyro_ = core::Vector3f::Zero();

    // bias correction
    core::Matrix3f dR_db_gyro_ = core::Matrix3f::Zero();
    core::Matrix3f dv_db_acc_  = core::Matrix3f::Zero();
    core::Matrix3f dv_db_gyro_ = core::Matrix3f::Zero();
    core::Matrix3f dp_db_acc_  = core::Matrix3f::Zero();
    core::Matrix3f dp_db_gyro_ = core::Matrix3f::Zero();

    float dT_; // total preintegration time
  };

  class ImuPreintegrator : public ImuPreintegratorBase {
  public:
    using CovType = Eigen::MatrixXf;

    using AType = Eigen::MatrixXf;
    using BType = Eigen::MatrixXf;

    void preintegrate(const ImuMeasurement& m, float dt) override;
    // allocates a new imu measurement
    void reset() override;
    // clang-format off
    inline const core::Matrix3f& delta_R() const override {return delta_R_;}
    inline const core::Vector3f& delta_p() const override {return delta_p_;}
    inline const core::Vector3f& delta_v() const override {return delta_v_;}
    inline const CovType& sigma() const {return sigma_;}
    // clang-format on

    // DEBUG function
    void getPrediction(const core::Isometry3f& Ti,
                       const core::Vector3f& vi,
                       core::Isometry3f& Tf,
                       core::Vector3f& vf);

    // protected:
    core::Vector3f delta_p_ = core::Vector3f::Zero();
    core::Matrix3f delta_R_ = core::Matrix3f::Identity();
    core::Vector3f delta_v_ = core::Vector3f::Zero();

    // allocate matrices for noise propagation
    AType A_ = AType::Identity(state_dim, state_dim);
    BType B_ = BType::Zero(state_dim, noise_dim);

    //
  };

} // namespace test_imu