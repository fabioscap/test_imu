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
    using CovType = core::MatrixN_<float, StateType::dim>;

    template <typename StateType>
    static void
    toMeanCov(const SigmaPoints<StateType>& spoints, StateType& mean, CovType<StateType>& cov) {
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
    static void toUnscented(const StateType& mean,
                            const CovType<StateType>& cov,
                            SigmaPoints<StateType>& spoints) {
      using TangentType = typename StateType::TangentType;

      const int state_dim = StateType::dim;
      const float lambda  = alpha_ * alpha_ * (state_dim + k_);

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

      Eigen::JacobiSVD<Eigen::MatrixXf> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
      double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
      std::cout << "det:  " << cov.determinant() << "\n";
      std::cout << "cond: " << cond << "\n";

      Eigen::MatrixXf U    = svd.matrixU();
      Eigen::MatrixXf S    = svd.singularValues();
      CovType<StateType> L = U * S.cwiseSqrt().asDiagonal();

      for (size_t i = 0; i < state_dim; ++i) {
        const TangentType Li = L.col(i);
        points[2 * i + 1]    = mean.boxplus(Li);
        points[2 * i + 2]    = mean.boxplus(-Li);
      }
    }

    // mysterious parameters
    static constexpr float beta_  = 2.0f;
    static constexpr float alpha_ = 1e-3;
    static constexpr float k_     = 0;

  private:
    UnscentedTransform();
  };

} // namespace test_imu