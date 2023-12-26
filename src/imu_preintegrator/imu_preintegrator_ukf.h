#pragma once

#include "common/manifold.h"
#include "common/unscented.h"

#include "imu_preintegrator_base.h"

namespace test_imu {

  /* The dynamic model of the IMU does not allow the input noises to be expressed additively */

  class ImuPreintegratorUKF : public ImuPreintegratorBase {
  public:
    using Vector3 = ImuPreintegratorBase::Vector3;
    using Matrix3 = ImuPreintegratorBase::Matrix3;

    static constexpr int state_dim = 15;
    static constexpr int input_dim = 6;

    using CovJType = core::MatrixN_<ImuPreintegratorBase::Scalar, state_dim + input_dim>;

    const core::Matrix3f delta_R() override;
    const core::Vector3f delta_p() override;
    const core::Vector3f delta_v() override;
    const core::MatrixX_<float> sigma() override;

    using DeltaManifold = ManifoldComp_<ManifoldSO3,   // dR
                                        Euclidean_<3>, // dv
                                        Euclidean_<3>, // dp
                                        Euclidean_<3>, // ba
                                        Euclidean_<3>, // bg
                                        Euclidean_<3>, // input acceleration
                                        Euclidean_<3>  // input gyroscope
                                        /* we add this as process noise
                                        Euclidean_<3>, // nba
                                        Euclidean_<3>  // nbg
                                        */
                                        >;

  protected:
    DeltaManifold delta_incr_;

    CovJType sigma_joint_                              = 1e-6 * CovJType::Identity();
    core::Vector_<Scalar, input_dim> input_noise_diag_ = core::Vector_<Scalar, input_dim>::Zero();
    core::Vector_<Scalar, input_dim> bias_noise_diag_  = core::Vector_<Scalar, input_dim>::Zero();

    // container for sigma points
    SigmaPoints<DeltaManifold> spoints;
    UnscentedTransform ut_;

    void preintegrate_(const ImuMeasurement& m, Scalar dt) override;
    void reset_() override;

    // dirty flag true iff delta_incr_ and sigma_joint_  are syncd with spoints
    bool is_updated_ = true;
  };

  class ImuPreintegratorUKFSlim : public ImuPreintegratorBase {
  public:
    using Vector3 = ImuPreintegratorBase::Vector3;
    using Matrix3 = ImuPreintegratorBase::Matrix3;

    static constexpr int state_dim = 9;
    static constexpr int input_dim = 6;

    using CovJType = core::MatrixN_<ImuPreintegratorBase::Scalar, state_dim + input_dim>;

    const core::Matrix3f delta_R() override;
    const core::Vector3f delta_p() override;
    const core::Vector3f delta_v() override;
    const core::MatrixX_<float> sigma() override;

    using DeltaManifold = ManifoldComp_<ManifoldSO3,   // dR
                                        Euclidean_<3>, // dv
                                        Euclidean_<3>, // dp
                                        Euclidean_<3>, // input acceleration
                                        Euclidean_<3>  // input gyroscope
                                        >;

  protected:
    DeltaManifold delta_incr_;

    CovJType sigma_joint_                              = 1e-6 * CovJType::Identity();
    core::Vector_<Scalar, input_dim> input_noise_diag_ = core::Vector_<Scalar, input_dim>::Zero();

    // container for sigma points
    SigmaPoints<DeltaManifold> spoints;
    UnscentedTransform ut_;

    void preintegrate_(const ImuMeasurement& m, Scalar dt) override;
    void reset_() override;

    // dirty flag true iff delta_incr_ and sigma_joint_  are syncd with spoints
    bool is_updated_ = true;

    // builds sigma_joint_ and geenrates sigma points
    void toUnscented_(float dt);
    // reconstruct mean and covariance and sets is_updated_ to true
    void toMeanCov_();
  };

} // namespace test_imu