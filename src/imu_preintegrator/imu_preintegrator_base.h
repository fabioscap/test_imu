#pragma once

#include "common/common.h"

namespace test_imu {

  // this holds the estimate jacobians wrt the biases.
  struct BiasJacobians {
    using Scalar = test_imu::Scalar; // float or double

    using Vector3 = core::Vector3_<Scalar>;
    using Matrix3 = core::Matrix3_<Scalar>;

    // bias correction
    Matrix3 dR_db_gyro = Matrix3::Zero();
    Matrix3 dv_db_acc  = Matrix3::Zero();
    Matrix3 dv_db_gyro = Matrix3::Zero();
    Matrix3 dp_db_acc  = Matrix3::Zero();
    Matrix3 dp_db_gyro = Matrix3::Zero();

    void update(const Matrix3& deltaR,
                const Matrix3& acc_skew,
                const Matrix3& dR,
                const Matrix3& Jr,
                Scalar dt);
  };

  class ImuPreintegratorBase {
  public:
    using Scalar = test_imu::Scalar; // float or double

    using Vector3 = core::Vector3_<Scalar>;
    using Matrix3 = core::Matrix3_<Scalar>;

    void preintegrate(const ImuMeasurement&, Scalar dt);
    void reset();

    virtual const core::Matrix3f delta_R()      = 0;
    virtual const core::Vector3f delta_v()      = 0;
    virtual const core::Vector3f delta_p()      = 0;
    virtual const core::MatrixX_<float> sigma() = 0;

    const void setNoiseGyroscope(const core::Vector3f& v);
    const void setNoiseAccelerometer(const core::Vector3f& v);
    const void setNoiseBiasGyroscope(const core::Vector3f& v);
    const void setNoiseBiasAccelerometer(const core::Vector3f& v);
    const void setBiasAcc(const core::Vector3f& v);
    const void setBiasGyro(const core::Vector3f& v);

    const core::Vector3f bias_gyro() const;
    const core::Vector3f bias_acc() const;

    std::vector<ImuMeasurement> measurements();
    const float dT() const;

    // forward step. does side effect on its arguments
    static void f(Matrix3& delta_R,
                  Vector3& delta_v,
                  Vector3& delta_p,
                  const Matrix3& dR,  // unbiased
                  const Vector3& acc, // unbiased
                  float dt);

    void getPrediction(const core::Isometry3f& Ti,
                       const core::Vector3f& vi,
                       core::Isometry3f& Tf,
                       core::Vector3f& vf);

    inline virtual const BiasJacobians& biasJacobians() const {
      return bias_J_;
    } // don't estimate biases by default

  protected:
    /* override these three functions */
    virtual void preintegrate_(const ImuMeasurement&, Scalar dt) = 0;
    virtual void reset_()                                        = 0;

    std::vector<ImuMeasurement> measurements_;
    Scalar dT_; // total preintegration time

    Vector3 acc_noise_       = Vector3::Zero();
    Vector3 gyro_noise_      = Vector3::Zero();
    Vector3 bias_acc_noise_  = Vector3::Zero();
    Vector3 bias_gyro_noise_ = Vector3::Zero();

    // nominal values for the bias
    Vector3 bias_acc_  = Vector3::Zero();
    Vector3 bias_gyro_ = Vector3::Zero();

    // some intermediate quantities that are used in two or more update functions
    Matrix3 dR_, Jr_;
    Vector3 acc_c_;
    float dt_;

    BiasJacobians bias_J_;
  };

} // namespace test_imu