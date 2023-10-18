#include "common/manifold.h"
#include "common/manifold_impl.cpp"
#include "common/unscented.h"
#include "common/unscented_impl.cpp"
#include <srrg_geometry/geometry3d.h>

using namespace test_imu;

int main() {
  using StateType = ManifoldSO3;
  using InputType = Euclidean_<3>;
  using JointType = ManifoldComp_<StateType, InputType>;

  using CovType = UnscentedTransform::CovType<JointType>;

  JointType s0;
  CovType cov = CovType::Identity();

  s0.get<0>().setData(Ry<Scalar>(0.3));
  s0.get<1>().setData(core::Vector3_<Scalar>(0.1, 0.2, 0.3));

  SigmaPoints<JointType> spoints;

  JointType rec;
  CovType cov_rec;

  UnscentedTransform::toUnscented(s0, cov, spoints);

  UnscentedTransform::toMeanCov(spoints, rec, cov_rec);

  std::cout << "s0\n";
  std::cout << s0.get<0>().data() << "\n" << s0.get<1>().data() << "\n";
  std::cout << cov << "\n";
  std::cout << "rec\n";
  std::cout << rec.get<0>().data() << "\n" << rec.get<1>().data() << "\n";
  std::cout << cov_rec << std::endl;
}