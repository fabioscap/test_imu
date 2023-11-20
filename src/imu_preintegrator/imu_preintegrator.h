#pragma once

#include "common/common.h"
#include "imu_preintegrator_base.h"

namespace test_imu {

  class ImuPreintegrator : public ImuPreintegratorBase {
  public:
    using CovType = ImuPreintegratorBase::CovType;

    using AType = Eigen::Matrix<ImuPreintegratorBase::Scalar, -1, -1>;
    using BType = Eigen::Matrix<ImuPreintegratorBase::Scalar, -1, -1>;

    using Matrix3 = ImuPreintegratorBase::Matrix3;
    using Vector3 = ImuPreintegratorBase::Vector3;

    void preintegrate(const ImuMeasurement& m, Scalar dt) override;
    // allocates a new imu measurement
    void reset() override;
    // clang-format off
    inline const core::Matrix3f delta_R() const override {return delta_R_.cast<float>();}
    inline const core::Vector3f delta_p() const override {return delta_p_.cast<float>();}
    inline const core::Vector3f delta_v() const override {return delta_v_.cast<float>();}
    inline const CovType& sigma() const {return sigma_;}
    // clang-format on

    // protected:
    Vector3 delta_p_ = Vector3::Zero();
    Matrix3 delta_R_ = Matrix3::Identity();
    Vector3 delta_v_ = Vector3::Zero();

    CovType sigma_ = 1e-10 * CovType::Identity(state_dim, state_dim);

    // allocate matrices for noise propagation
    AType A_ = AType::Identity(state_dim, state_dim);
    BType B_ = BType::Zero(state_dim, noise_dim);

    //
  };

  // Imu Preintegrator version with no bias estimation and bias correction
  class ImuPreintegratorSlim : public ImuPreintegratorBase {
  public:
    using CovType = ImuPreintegratorBase::CovType;

    using AType = Eigen::Matrix<ImuPreintegratorBase::Scalar, -1, -1>;
    using BType = Eigen::Matrix<ImuPreintegratorBase::Scalar, -1, -1>;

    using Matrix3 = ImuPreintegratorBase::Matrix3;
    using Vector3 = ImuPreintegratorBase::Vector3;

    static constexpr int state_dim = 9;
    static constexpr int noise_dim = 6;

    void preintegrate(const ImuMeasurement& m, Scalar dt) override;
    // allocates a new imu measurement
    void reset() override;
    // clang-format off
    inline const core::Matrix3f delta_R() const override {return delta_R_.cast<float>();}
    inline const core::Vector3f delta_p() const override {return delta_p_.cast<float>();}
    inline const core::Vector3f delta_v() const override {return delta_v_.cast<float>();}
    inline const CovType& sigma() const {return sigma_;}
    // clang-format on

  protected:
    Vector3 delta_p_ = Vector3::Zero();
    Matrix3 delta_R_ = Matrix3::Identity();
    Vector3 delta_v_ = Vector3::Zero();

    CovType sigma_ = 1e-10 * CovType::Identity(state_dim, state_dim);

    // allocate matrices for noise propagation
    AType A_ = AType::Identity(state_dim, state_dim);
    BType B_ = BType::Zero(state_dim, noise_dim);
  };

} // namespace test_imu