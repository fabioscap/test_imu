#include "common/manifold.h"
#include "common/ukf.h"
#include <srrg_geometry/geometry3d.h>

using namespace test_imu;

int main() {
  using StateType = ManifoldSO3;
  using CovType   = UnscentedTransform::CovType<StateType>;

  StateType mean = StateType();
  CovType cov    = CovType::Identity();

  mean.setData(Ry(0.2f));
  cov(4, 4) = 0.5;

  SigmaPoints<StateType> spoints;
  UnscentedTransform::toUnscented(mean, cov, spoints);

  StateType mean_rec;
  CovType cov_rec;

  UnscentedTransform::toMeanCov(spoints, mean_rec, cov_rec);

  std::cout << "mean: \n"
            << mean.data() << "\n"
            << "cov: \n"
            << cov << "\n";
  std::cout << "mean rec: \n"
            << mean_rec.data() << "\n"
            << "cov rec: \n"
            << cov_rec << "\n";
}