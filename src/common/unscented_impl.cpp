#include "unscented.h"

namespace test_imu {

  // utility functions
  template <typename StateType>
  void UnscentedTransform::compute_weights_scaling(SigmaPoints<StateType>& spoints,
                                                   Scalar& cov_scaling) const {
    constexpr int state_dim = StateType::dim;

    Scalar& wm0 = spoints.wm0;
    Scalar& wc0 = spoints.wc0;
    Scalar& wmi = spoints.wmi;
    Scalar& wci = spoints.wci;

    if (weight_scheme_ == WeightScheme::UKF) {
      const Scalar lambda = alpha_ * alpha_ * (state_dim + k_) - state_dim;

      // weights 0
      wm0 = lambda / (state_dim + lambda);
      wc0 = wm0 + (1 - alpha_ * alpha_ + beta_);
      // weights i
      wmi = 0.5 / Scalar(state_dim + lambda);
      wci = 0.5 / Scalar(state_dim + lambda);

      cov_scaling = (state_dim + lambda);
    } else if (weight_scheme_ == WeightScheme::HB) {
      /* hertzberg weights*/
      wm0         = 1.0 / (2 * state_dim + 1);
      wmi         = wm0;
      wc0         = 0.5;
      wci         = wc0;
      cov_scaling = 1.0;
    } else {
      throw std::runtime_error(
        "UnscentedTransform::compute_weights_scaling| unknown weight scheme");
    }
  }

  template <typename StateType>
  void UnscentedTransform::compute_sigma_points(const StateType& mean,
                                                const CovType<StateType>& L_scaled,
                                                SigmaPoints<StateType>& spoints) const {
    constexpr int state_dim = StateType::dim;
    using TangentType       = typename StateType::TangentType;
    auto& points            = spoints.points;
    points[0]               = mean;
    for (size_t i = 0; i < state_dim; ++i) {
      const TangentType Li = L_scaled.col(i);
      points[2 * i + 1]    = mean.boxplus(Li);
      points[2 * i + 2]    = mean.boxplus(-Li);
    }
  }

  template <typename StateType>
  void UnscentedTransform::mean_from_sigma_points(const SigmaPoints<StateType>& spoints,
                                                  StateType& mean,
                                                  const size_t n_iters) const {
    using TangentType  = typename StateType::TangentType;
    const auto& points = spoints.points;
    const Scalar wm0   = spoints.wm0;
    const Scalar wmi   = spoints.wmi;

    // compute chart point
    mean = points[0];

    for (size_t i = 0; i < n_iters; ++i) {
      // compute tangent mean displacement
      TangentType displacement = wm0 * mean.boxminus(points.at(0));
      for (size_t i = 1; i < points.size(); ++i) {
        const StateType& Xi = points.at(i);
        displacement += wmi * mean.boxminus(Xi);
      }
      // add it to the mean
      mean = mean.boxplus(displacement);
    }
  }

  template <typename StateType>
  void UnscentedTransform::toMeanCov(const SigmaPoints<StateType>& spoints,
                                     StateType& mean,
                                     CovType<StateType>& cov) const {
    using TangentType  = typename StateType::TangentType;
    const auto& points = spoints.points;
    const Scalar wc0   = spoints.wc0;
    const Scalar wci   = spoints.wci;
    mean_from_sigma_points(spoints, mean, 0);

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
                                       SigmaPoints<StateType>& spoints) const {
    constexpr int state_dim = StateType::dim;
    auto& points            = spoints.points;

    // points.clear();
    points.resize(2 * state_dim + 1);

    Scalar cov_scaling;

    compute_weights_scaling(spoints, cov_scaling);
    // Perform Cholesky decomposition
    // TODO use selfAdjointView where it is possible
    Eigen::LLT<CovType<StateType>> llt(cov_scaling *
                                       (cov + CovType<StateType>::Identity() * cov_regularizer_));

    if (!llt.info() == Eigen::Success)
      throw std::runtime_error("UnscentedTransform::toUnscented| Cholesky decomposition failed");

    const CovType<StateType>& L = llt.matrixL();
    compute_sigma_points(mean, L, spoints);
  }

  template <typename StateType>
  void UnscentedTransform::toMeanSqrtCov(const SigmaPoints<StateType>& spoints,
                                         StateType& mean,
                                         CovType<StateType>& cov_sqrt) const {
    constexpr int state_dim = StateType::dim;
    const auto& points      = spoints.points;
    const Scalar wc0        = spoints.wc0;
    const Scalar wci        = spoints.wci;

    mean_from_sigma_points(spoints, mean, 0);

    // QR

    // maybe I can preallocate this stupid matrix
    Eigen::Matrix<Scalar, state_dim, 2 * state_dim> C;

    Scalar wci_sqrt = std::sqrt(wci);
    for (size_t i = 1; i < points.size(); ++i) {
      C.block(0, i - 1, state_dim, 1) = wci_sqrt * mean.boxminus(points.at(i));
    }

    cov_sqrt = C.transpose()
                 .householderQr()
                 .matrixQR()
                 .topLeftCorner(state_dim, state_dim)
                 .template triangularView<Eigen::Upper>();

    // there is no need to redo llt here, I should implement a choleskyUpdate function that
    // operates directly on cov_sqrt
    Eigen::internal::llt_inplace<Scalar, Eigen::Upper>::rankUpdate(
      cov_sqrt, mean.boxminus(points.at(0)), wc0);
    cov_sqrt.transposeInPlace();
  }

  template <typename StateType>
  void UnscentedTransform::toSqrtUnscented(const StateType& mean,
                                           const CovType<StateType>& cov_sqrt,
                                           SigmaPoints<StateType>& spoints) const {
    constexpr int state_dim = StateType::dim;

    auto& points = spoints.points;
    points.resize(2 * state_dim + 1);

    Scalar cov_scaling;

    compute_weights_scaling(spoints, cov_scaling);
    compute_sigma_points(mean, static_cast<CovType<StateType>>(cov_scaling * cov_sqrt), spoints);
  }
} // namespace test_imu