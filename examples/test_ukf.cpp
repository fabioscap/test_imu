#include "common/manifold.h"
#include "common/ukf.h"
#include <srrg_geometry/geometry3d.h>

using namespace test_imu;

int main() {
  using UKFType = UKF_<ManifoldComp_<ManifoldSO3, Euclidean_<3>>>;

  UKFType ukf;

  UKFType::StateType mean = UKFType::StateType();

  mean.get<0>().setData(Ry(0.2f));
  mean.get<1>().setData(core::Vector3f(1.0f, 2.0f, 3.0f));

  UKFType::CovType cov = UKFType::CovType::Identity();

  cov(4, 4) = 0.5;

  ukf.toUnscented(mean, cov);

  UKFType::StateType mean_rec;
  UKFType::CovType cov_rec;

  ukf.toMeanCov(mean_rec, cov_rec);

  std::cout << "mean: \n"
            << mean.get<0>().data() << "\n"
            << mean.get<1>().data() << "\n---\n"
            << "cov: \n"
            << cov << "\n";
  std::cout << "mean rec: \n"
            << mean_rec.get<0>().data() << "\n"
            << mean_rec.get<1>().data() << "\n---\n"
            << "cov rec: \n"
            << cov_rec << "\n";
}