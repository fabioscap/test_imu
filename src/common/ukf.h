#include "common.h"

#include <Eigen/Dense>

namespace test_imu {

  // Only what is needed for prediction. srrg2_solver does the update.
  template <typename ManifoldType_>
  class UKF_ {
  public:
    static constexpr int state_dim = ManifoldType_::dim;
    static constexpr int N         = 2 * state_dim + 1;

    using StateType   = ManifoldType_;
    using TangentType = typename ManifoldType_::TangentType;
    using CovType     = core::MatrixN_<float, state_dim>;

    void toMeanCov(StateType& mean, CovType& cov) const {
      // compute chart point
      const StateType& sigma_0 = points_[0];

      // compute tangent mean displacement
      TangentType displacement = TangentType::Zero();
      for (size_t i = 1; i < points_.size(); ++i) {
        const StateType& Xi = points_.at(i);
        displacement += wi_ * sigma_0.boxminus(Xi);
      }

      // add it to the mean
      mean = sigma_0.boxplus(displacement);

      // covariance lives in tangent space
      cov.setZero();
      TangentType err = mean.boxminus(points_.at(0));
      cov += wc0_ * err * err.transpose();

      for (size_t i = 1; i < points_.size(); ++i) {
        TangentType err = mean.boxminus(points_.at(i));
        cov += wi_ * err * err.transpose();
      }
    }
    void toUnscented(const StateType& mean, const CovType& cov) {
      points_[0] = mean;

      // mean weight
      wm0_ = lambda_ / (state_dim + lambda_);
      // covariance weight
      wc0_ = wm0_ + (1 - alpha_ * alpha_ + beta_);
      // weight i
      wi_ = 0.5 / (state_dim + lambda_);

      // Perform Cholesky decomposition
      Eigen::LLT<Eigen::MatrixXf> llt((state_dim + lambda_) * cov);

      if (!llt.info() == Eigen::Success)
        throw std::runtime_error("UKF::toUnscented| Cholesky decomposition failed");

      CovType L = llt.matrixL(); // should be symmetric
      for (size_t i = 0; i < state_dim; ++i) {
        const TangentType Li = L.col(i);
        points_[2 * i + 1]   = mean.boxplus(Li);
        points_[2 * i + 2]   = mean.boxplus(-Li);
      }
    }

    // protected:
    std::vector<StateType> points_ = std::vector<StateType>(N);
    float wm0_, wc0_, wi_;

    // mysterious parameters
    float beta_  = 2.0f;
    float alpha_ = 1e-3;
    float k_     = 0;

    float lambda_ = alpha_ * alpha_ * (state_dim + k_);
  };

} // namespace test_imu