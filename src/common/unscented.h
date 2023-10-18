#pragma once

#include "common.h"

#include <Eigen/Dense>

namespace test_imu {

  template <typename StateType>
  struct SigmaPoints {
    std::vector<StateType> points;
    float wm0, wc0, wi;
  };

  class UnscentedTransform {
  public:
    template <typename StateType>
    using CovType = core::MatrixN_<Scalar, StateType::dim>;

    template <typename StateType>
    static void
    toMeanCov(const SigmaPoints<StateType>& spoints, StateType& mean, CovType<StateType>& cov);

    template <typename StateType>
    static void toUnscented(const StateType& mean,
                            const CovType<StateType>& cov,
                            SigmaPoints<StateType>& spoints);

    // mysterious parameters
    static constexpr float beta_  = 2.0f;
    static constexpr float alpha_ = 2e-2;
    static constexpr float k_     = 0;

  private:
    UnscentedTransform();
  };

} // namespace test_imu