#include "common/manifold.h"
#include "common/ukf.h"
#include <srrg_geometry/geometry3d.h>

using namespace test_imu;

int main() {
  using UKFType = UKF_<ManifoldSO3>;

  UKFType ukf;

  UKFType::StateType mean = UKFType::StateType(Ry(0.2f));

  UKFType::CovType cov = UKFType::CovType::Identity();

  ukf.toUnscented(mean, cov);

  UKFType::StateType mean_rec;
  UKFType::CovType cov_rec;

  ukf.toMeanCov(mean_rec, cov_rec);

  std::cout << "mean: \n" << mean.data() << "\ncov: \n" << cov << "\n";
  std::cout << "mean rec: \n" << mean_rec.data() << "\ncov rec: \n" << cov_rec << "\n";
}