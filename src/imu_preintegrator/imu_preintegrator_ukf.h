#pragma once

#include "common/common.h"
#include "common/manifold.h"
#include "common/unscented.h"

#include "imu_preintegrator.h"

namespace test_imu {

  class ImuPreintegratorUKF : public ImuPreintegratorBase {
  public:
    using CovType      = ImuPreintegratorBase::CovType;
    using CovNoiseType = ImuPreintegratorBase::CovNoiseType;

    using Vector3 = ImuPreintegratorBase::Vector3;
    using Matrix3 = ImuPreintegratorBase::Matrix3;

    static constexpr int state_dim = 15;
    static constexpr int noise_dim = 12;

    using CovJType = core::MatrixN_<ImuPreintegratorBase::Scalar, state_dim + noise_dim>;

    void preintegrate(const ImuMeasurement& m, Scalar dt) override;
    // allocates a new imu measurement
    void reset() override;
    // clang-format off
    const core::Matrix3f delta_R() const override;
    const core::Vector3f delta_p() const override;
    const core::Vector3f delta_v() const override;
    inline const CovType& sigma() const override {
      // sigma_ = sigma_joint_.block<state_dim, state_dim>(0, 0);
      return sigma_;}
    // clang-format on

    using DeltaManifold = ManifoldComp_<ManifoldSO3,   // dR
                                        Euclidean_<3>, // dv
                                        Euclidean_<3>, // dp
                                        Euclidean_<3>, // ba
                                        Euclidean_<3>, // bg
                                        Euclidean_<3>, // na
                                        Euclidean_<3>, // ng
                                        Euclidean_<3>, // nba
                                        Euclidean_<3>  // nbg
                                        >;
    // DEBUG function
    void getPrediction(const core::Isometry3f& Ti,
                       const core::Vector3f& vi,
                       core::Isometry3f& Tf,
                       core::Vector3f& vf);

    // protected:
    DeltaManifold delta_incr_;

    // UKF: preallocate sigma joint input and state
    CovJType sigma_joint_ = 1e-6 * CovJType::Identity();

    // container for sigma points
    SigmaPoints<DeltaManifold> spoints;
  };

} // namespace test_imu