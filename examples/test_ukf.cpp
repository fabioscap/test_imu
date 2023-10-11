#include "common/ukf.h"
#include <srrg_geometry/geometry3d.h>

using namespace test_imu;

using UKFType = UKF<core::Matrix3f, 3>;

core::Matrix3f boxplusSO3(const core::Matrix3f& from, const core::Vector3f& dsp) {
  core::Matrix3f R = from * core::geometry3d::expMapSO3(dsp);
  core::fixRotation(R);
  return R;
}
core::Vector3f boxminusSO3(const core::Matrix3f& from, const core::Matrix3f& to) {
  core::Matrix3f R = from.transpose() * to;
  return core::geometry3d::logMapSO3(R.eval());
}
core::Vector3f boxplusVector(const core::Vector3f& from, const core::Vector3f& dsp) {
  return from + dsp;
}
core::Vector3f boxminusVector(const core::Vector3f& from, const core::Vector3f& to) {
  return from - to;
}

int main() {
  // to unscented and back

  /*   UKFType ukf(boxplusSO3, boxminusSO3);

    UKFType::StateType mean = Ry(0.3f);

    UKFType::CovType cov = UKFType::CovType::Identity();
    ukf.toUnscented(mean, cov);

    UKFType::StateType mean_rec;
    UKFType::CovType cov_rec;
    ukf.toMeanCov(mean_rec, cov_rec);

    std::cout << "original: \n" << mean << "\n---\n" << cov << "\n" << std::endl;
    std::cout << "reconstructed: \n" << mean_rec << "\n---\n" << cov_rec << std::endl; */
}