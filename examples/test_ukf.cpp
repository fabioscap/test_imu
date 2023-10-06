#include "common/ukf.h"
#include <srrg_geometry/geometry3d.h>

using namespace test_imu;

using UKFType = UKF<core::Matrix3f, 3>;

namespace test_imu {
  template <>
  UKFType::StateType boxplus(const UKFType::StateType& from, const UKFType::TangentType& dsp) {
    core::Matrix3f R = from * core::geometry3d::expMapSO3(dsp);
    core::fixRotation(R);
    return R;
  }
  template <>
  UKFType::TangentType boxminus(const UKFType::StateType& from, const UKFType::StateType& to) {
    core::Matrix3f R = from.transpose() * to;
    return core::geometry3d::logMapSO3(R.eval());
  }
} // namespace test_imu

int main() {
  // to unscented and back

  UKFType ukf;

  UKFType::StateType mean = UKFType::StateType::Identity();

  UKFType::CovType cov = UKFType::CovType::Identity();
  cov(0, 0)            = 10;

  ukf.toUnscented(mean, cov);

  UKFType::StateType mean_rec = UKFType::StateType::Identity();
  UKFType::CovType cov_rec    = UKFType::CovType::Identity();
  ukf.toMeanCov(mean, cov);

  std::cout << "sigma points:\n";
  for (size_t i = 0; i < ukf.points_.size(); ++i) {
    std::cout << "-" << i << "-\n" << ukf.points_.at(i) << "\n";
  }

  std::cout << "mean: \n";
  std::cout << mean_rec << "\n";
  std::cout << "covariance: \n";
  std::cout << cov_rec << "\n";

  // pass through a function
}