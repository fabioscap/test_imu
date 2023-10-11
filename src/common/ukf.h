#include "common.h"
#include <Eigen/Dense>

namespace test_imu {

  template <int dim>
  using TangentType_ = core::Vector_<float, dim>;

  // Only what is needed for prediction. srrg2_solver does the update.

  template <typename StateType_, int state_dim_>
  class UKF {
  public:
    static constexpr int state_dim = state_dim_;
    static constexpr int N         = 2 * state_dim + 1;

    using StateType   = StateType_;
    using TangentType = core::Vector_<float, state_dim>;
    using CovType     = core::MatrixN_<float, state_dim>;

    void toMeanCov(StateType_& mean, CovType& cov) const {
      // compute chart point
      const StateType_& sigma_0 = points_[0];

      // compute tangent mean displacement
      TangentType displacement = TangentType::Zero();
      for (size_t i = 1; i < points_.size(); ++i) {
        const StateType_& Xi = points_.at(i);
        displacement -= wi_ * boxminus(sigma_0, Xi);
      }

      // add it to the mean
      mean = boxplus(sigma_0, displacement);

      // covariance lives in tangent space
      cov.setZero();
      TangentType err = boxminus(mean, points_.at(0));
      cov += wc0_ * err * err.transpose();

      for (size_t i = 1; i < points_.size(); ++i) {
        TangentType err = boxminus(mean, points_.at(i));
        cov += wi_ * err * err.transpose();
      }
    }
    void toUnscented(const StateType_& mean, const CovType& cov) {
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
        points_[2 * i + 1]   = boxplus(mean, Li);
        points_[2 * i + 2]   = boxplus(mean, -Li);
      }
    }

    // protected:
    std::vector<StateType_> points_ = std::vector<StateType_>(N);
    float wm0_, wc0_, wi_;

    // mysterious parameters
    float beta_  = 2.0f;
    float alpha_ = 1e-3;
    float k_     = 0;

    float lambda_ = alpha_ * alpha_ * (state_dim + k_);
  };

  template <typename DataType_, int dim_>
  struct ManifoldBase {
    static constexpr int dim = dim_;
    using DataType           = DataType_;
    using TangentType        = TangentType_<dim>;
  };

} // namespace test_imu