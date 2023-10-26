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
    void toMeanCov(const SigmaPoints<StateType>& spoints, StateType& mean, CovType<StateType>& cov);

    template <typename StateType>
    void toUnscented(const StateType& mean,
                     const CovType<StateType>& cov,
                     SigmaPoints<StateType>& spoints);

    // protected:
    // mysterious parameters
    float beta_  = 2.0f;
    float alpha_ = 2e-2;
    float k_     = 0;

    WeightScheme weight_scheme_ = WeightScheme::UKF;

    float cov_regularizer_ = 1e-10;
    /*   private:
        UnscentedTransform(); */
  };

} // namespace test_imu