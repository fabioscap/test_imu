#include "unscented.h"

namespace test_imu {

  template <typename StateType>
  void UnscentedTransform::toMeanCov(const SigmaPoints<StateType>& spoints,
                                     StateType& mean,
                                     CovType<StateType>& cov) {
    using TangentType = typename StateType::TangentType;

    const auto& points = spoints.points;
    const float wc0    = spoints.wc0;
    const float wmi    = spoints.wmi;
    const float wci    = spoints.wci;

    // compute chart point
    mean        = points[0];
    int n_iters = 100;
    for (int i = 0; i < n_iters; ++i) { // TODO do this iteration until mean does not change
      // compute tangent mean displacement
      TangentType displacement = spoints.wm0 * mean.boxminus(points.at(0));
      for (size_t i = 1; i < points.size(); ++i) {
        const StateType& Xi = points.at(i);
        displacement += wmi * mean.boxminus(Xi);
      }
      // add it to the mean
      mean = mean.boxplus(displacement);
    }

    // covariance lives in tangent space
    cov.setZero();
    TangentType err = mean.boxminus(points.at(0));
    cov += wc0 * err * err.transpose();

    for (size_t i = 1; i < points.size(); ++i) {
      err = mean.boxminus(points.at(i));
      cov += wci * err * err.transpose();
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
    float& wmi   = spoints.wmi;
    float& wci   = spoints.wci;

    points.resize(2 * state_dim + 1);

    points[0] = mean;

    // mean weight
    /* wm0 = lambda / (state_dim + lambda);
    // covariance weight
    wc0 = wm0 + (1 - alpha_ * alpha_ + beta_);
    // weight i
    wmi = 0.5 / (state_dim + lambda);

    wci = 0.5 / (state_dim + lambda); */

    /* hertzberg weights*/
    wm0 = 1 / (2 * state_dim + 1);
    wmi = wm0;
    wc0 = 0.5;
    wci = 0.5;

    /*       Eigen::JacobiSVD<Eigen::MatrixXf> svd(cov, Eigen::ComputeFullU |
       Eigen::ComputeFullV); Eigen::MatrixXf U = svd.matrixU(); Eigen::MatrixXf S =
       svd.singularValues(); CovType<StateType> L = U * std::sqrt((state_dim + lambda)) *
       S.cwiseSqrt().asDiagonal() * U.transpose(); */

    float cov_regularizer = 1e-10;

    // Perform Cholesky decomposition
    // Eigen::LLT<CovType<StateType>> llt((state_dim + lambda) *
    //                                   (cov + CovType<StateType>::Identity() * cov_regularizer));

    // hertzberg weights
    Eigen::LLT<CovType<StateType>> llt((cov + CovType<StateType>::Identity() * cov_regularizer));

    // Eigen::LLT<CovType<StateType>> llt((cov + CovType<StateType>::Identity() * cov_regularizer));

    if (!llt.info() == Eigen::Success)
      throw std::runtime_error("UnscentedTransform::toUnscented| Cholesky decomposition failed");

    CovType<StateType> L = llt.matrixL();

    for (size_t i = 0; i < state_dim; ++i) {
      const TangentType Li = L.col(i);
      points[2 * i + 1]    = mean.boxplus(Li);
      points[2 * i + 2]    = mean.boxplus(-Li);
    }
  }
} // namespace test_imu