#include "unscented.h"

namespace test_imu {

  template <typename StateType>
  void UnscentedTransform::toMeanCov(const SigmaPoints<StateType>& spoints,
                                     StateType& mean,
                                     CovType<StateType>& cov) {
    using TangentType = typename StateType::TangentType;

    const auto& points = spoints.points;
    const float wc0    = spoints.wc0;
    const float wi     = spoints.wi;

    // compute chart point
    const StateType& sigma_0 = points[0];

    // compute tangent mean displacement
    TangentType displacement = TangentType::Zero();
    for (size_t i = 1; i < points.size(); ++i) {
      const StateType& Xi = points.at(i);
      displacement += wi * sigma_0.boxminus(Xi);
    }

    // add it to the mean
    mean = sigma_0.boxplus(displacement);

    // covariance lives in tangent space
    cov.setZero();
    TangentType err = mean.boxminus(sigma_0);
    cov += wc0 * err * err.transpose();

    for (size_t i = 1; i < points.size(); ++i) {
      TangentType err = mean.boxminus(points.at(i));
      cov += wi * err * err.transpose();
    }
  }

  template <typename StateType>
  void UnscentedTransform::toUnscented(const StateType& mean,
                                       const CovType<StateType>& cov,
                                       SigmaPoints<StateType>& spoints) {
    using TangentType = typename StateType::TangentType;

    const int state_dim = StateType::dim;
    const float lambda  = alpha_ * alpha_ * (state_dim + k_) - state_dim;

    auto& points = spoints.points;
    float& wm0   = spoints.wm0;
    float& wc0   = spoints.wc0;
    float& wi    = spoints.wi;

    points.resize(2 * state_dim + 1);

    points[0] = mean;

    // mean weight
    wm0 = lambda / (state_dim + lambda);
    // covariance weight
    wc0 = wm0 + (1 - alpha_ * alpha_ + beta_);
    // weight i
    wi = 0.5 / (state_dim + lambda);

    /*       Eigen::JacobiSVD<Eigen::MatrixXf> svd(cov, Eigen::ComputeFullU |
       Eigen::ComputeFullV); Eigen::MatrixXf U = svd.matrixU(); Eigen::MatrixXf S =
       svd.singularValues(); CovType<StateType> L = U * std::sqrt((state_dim + lambda)) *
       S.cwiseSqrt().asDiagonal() * U.transpose(); */

    CovType<StateType> L = (cov).llt().matrixL();

    for (size_t i = 0; i < state_dim; ++i) {
      const TangentType Li = L.col(i);
      points[2 * i + 1]    = mean.boxplus(std::sqrt(state_dim + lambda) * Li);
      points[2 * i + 2]    = mean.boxplus(-std::sqrt(state_dim + lambda) * Li);
    }
  }
} // namespace test_imu