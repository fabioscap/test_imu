#pragma once

#include "common.h"

#include <Eigen/Dense>

namespace test_imu {

  enum class WeightScheme { UKF, HB };

  template <typename StateType>
  struct SigmaPoints {
    std::vector<StateType> points;
    float wm0, wc0, wmi, wci;
  };

  class UnscentedTransform {
  public:
    template <typename StateType>
    using CovType = core::MatrixN_<Scalar, StateType::dim>;

    template <typename StateType>
    void toMeanCov(const SigmaPoints<StateType>& spoints,
                   StateType& mean,
                   CovType<StateType>& cov) const;

    template <typename StateType>
    void toUnscented(const StateType& mean,
                     const CovType<StateType>& cov,
                     SigmaPoints<StateType>& spoints) const;

    // for square root UKF (van der Merwe), keep covariance in square root form
    template <typename StateType>
    void toMeanSqrtCov(const SigmaPoints<StateType>& spoints,
                       const CovType<StateType>& process_noise_cov_sqrt,
                       StateType& mean,
                       CovType<StateType>& cov_sqrt) const;
    template <typename StateType>
    void toSqrtUnscented(const StateType& mean,
                         const CovType<StateType>& cov_sqrt,
                         SigmaPoints<StateType>& spoints) const;

    // protected:
    // mysterious parameters
    float beta_  = 2.0f;
    float alpha_ = 1e-3;
    float k_     = 0;

    WeightScheme weight_scheme_ = WeightScheme::UKF;

    float cov_regularizer_ = 1e-10;

  private:
    template <typename StateType>
    void compute_weights_scaling(SigmaPoints<StateType>& spoints, float& cov_scaling) const;
    template <typename StateType>
    void compute_sigma_points(const StateType& mean,
                              const CovType<StateType>& L_scaled,
                              SigmaPoints<StateType>& spoints) const;
    template <typename StateType>
    void mean_from_sigma_points(const SigmaPoints<StateType>& spoints,
                                StateType& mean,
                                const size_t n_iters = 10) const;
  };

} // namespace test_imu