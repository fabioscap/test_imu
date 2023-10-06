#include "common.h"
#include <Eigen/Dense>

namespace test_imu {

  template <typename StateType, typename TangentType>
  inline TangentType boxminus(const StateType& from, const StateType& to) {
    throw std::runtime_error("boxminus not implemented");
  }
  template <typename StateType, typename TangentType>
  inline StateType boxplus(const StateType& from, const TangentType& disp) {
    throw std::runtime_error("boxplus not implemented");
  }

  template <typename StateType_, int state_dim_>
  class UKF {
  public:
    static constexpr int state_dim = state_dim_;
    static constexpr int N         = 2 * state_dim + 1;

    using StateType   = StateType_;
    using TangentType = core::Vector_<float, state_dim>;
    using CovType     = core::MatrixN_<float, state_dim>;

    void toMeanCov(StateType& mean, CovType& cov) const {
      // compute chart point
      const StateType& sigma_0 = points_[0];

      // compute tangent mean displacement
      TangentType displacement = TangentType::Zero();
      for (size_t i = 1; i < points_.size(); ++i)
        displacement -= w_i_ * boxminus<StateType, TangentType>(sigma_0, points_.at(i));

      // add it to the mean
      mean = boxplus(sigma_0, displacement);

      // covariance lives in tangent space
      for (size_t i = 0; i < points_.size(); ++i) {
        const TangentType err = boxminus<StateType, TangentType>(mean, points_.at(i));
        cov += err * err.transpose();
      }
    }
    void toUnscented(const StateType& mean, const CovType& cov) {
      points_[0] = mean;

      // mean weight
      w_m_0 = lambda_ / (state_dim + lambda_);
      // covariance weight
      w_c_0 = w_m_0 + (1 - alpha_ * alpha_ + beta_);
      // weight i
      w_i_ = 0.5 / (state_dim + lambda_);

      // perform cholesky decomposition
      // Perform Cholesky decomposition
      Eigen::LLT<Eigen::MatrixXf> llt((state_dim + lambda_) * cov);

      if (!llt.info() == Eigen::Success)
        throw std::runtime_error("UKF::toUnscented| cholesky decomposition failed");

      CovType L = llt.matrixL(); // should be symmetric
      for (size_t i = 0; i < state_dim; ++i) {
        const TangentType Li = L.col(i);
        points_[2 * i + 1]   = boxplus(mean, Li);
        points_[2 * i + 2]   = boxplus(mean, Li);
      }
    }

    // protected:
    std::vector<StateType> points_ = std::vector<StateType>(N);
    float w_m_0, w_c_0, w_i_;

    // mysterious parameters
    float beta_  = 2.0f;
    float alpha_ = 1e-3;
    float k_     = 0;

    float lambda_ = alpha_ * alpha_ * state_dim;
  };

} // namespace test_imu