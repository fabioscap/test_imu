#pragma once

#include "common/common.h"
#include "imu_preintegrator_base.h"

namespace test_imu {

  class ImuPreintegrator : public ImuPreintegratorBase {
  public:
    using Matrix3 = ImuPreintegratorBase::Matrix3;
    using Vector3 = ImuPreintegratorBase::Vector3;

    static constexpr int state_dim = 15;
    static constexpr int noise_dim = 12;

    using CovType = core::MatrixN_<Scalar, state_dim>;

    using AType = Eigen::Matrix<ImuPreintegratorBase::Scalar, state_dim, state_dim>;
    using BType = Eigen::Matrix<ImuPreintegratorBase::Scalar, state_dim, noise_dim>;

    // clang-format off
    inline const core::Matrix3f delta_R() override {return delta_R_.cast<float>();}
    inline const core::Vector3f delta_p() override {return delta_p_.cast<float>();}
    inline const core::Vector3f delta_v() override {return delta_v_.cast<float>();}
    inline const core::MatrixX_<float> sigma() override {return sigma_.cast<float>();}
    // clang-format on

  protected:
    Vector3 delta_p_ = Vector3::Zero();
    Matrix3 delta_R_ = Matrix3::Identity();
    Vector3 delta_v_ = Vector3::Zero();

    CovType sigma_ = 1e-10 * CovType::Identity(state_dim, state_dim);
    core::Vector_<Scalar, noise_dim> sigma_noise_diag_ = core::Vector_<Scalar, noise_dim>::Zero();

    // allocate matrices for noise propagation
    AType A_ = AType::Identity(state_dim, state_dim);
    BType B_ = BType::Zero(state_dim, noise_dim);
    //

    void preintegrate_(const ImuMeasurement& m, Scalar dt) override;
    void reset_() override;
  };

  // Imu Preintegrator version with no bias estimation and bias correction
  class ImuPreintegratorSlim : public ImuPreintegratorBase {
  public:
    using Matrix3 = ImuPreintegratorBase::Matrix3;
    using Vector3 = ImuPreintegratorBase::Vector3;

    static constexpr int state_dim = 9;
    static constexpr int noise_dim = 6;

    using CovType = core::MatrixN_<Scalar, state_dim>;

    using AType = Eigen::Matrix<ImuPreintegratorBase::Scalar, state_dim, state_dim>;
    using BType = Eigen::Matrix<ImuPreintegratorBase::Scalar, state_dim, noise_dim>;

    // clang-format off
    inline const core::Matrix3f delta_R() override {return delta_R_.cast<float>();}
    inline const core::Vector3f delta_p() override {return delta_p_.cast<float>();}
    inline const core::Vector3f delta_v() override {return delta_v_.cast<float>();}
    inline const core::MatrixX_<float> sigma() override {return sigma_.cast<float>();}
    // clang-format on

  protected:
    Vector3 delta_p_ = Vector3::Zero();
    Matrix3 delta_R_ = Matrix3::Identity();
    Vector3 delta_v_ = Vector3::Zero();

    CovType sigma_ = 1e-10 * CovType::Identity(state_dim, state_dim);
    core::Vector_<Scalar, noise_dim> sigma_noise_diag_ = core::Vector_<Scalar, noise_dim>::Zero();

    // allocate matrices for noise propagation
    AType A_ = AType::Identity(state_dim, state_dim);
    BType B_ = BType::Zero(state_dim, noise_dim);
    //

    void preintegrate_(const ImuMeasurement& m, Scalar dt) override;
    void reset_() override;
  };

} // namespace test_imu