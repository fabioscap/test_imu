#include "common/manifold.h"
#include "common/manifold_impl.cpp"
#include "common/unscented.h"
#include "common/unscented_impl.cpp"
#include <srrg_geometry/geometry3d.h>

#include <string.h>
using namespace test_imu;

using StateType = ManifoldSO2;
using InputType = Euclidean_<1>;
using JointType = ManifoldComp_<StateType, InputType>;

using CovType = UnscentedTransform::CovType<JointType>;

JointType f(const JointType& x) {
  JointType x_next;

  core::Matrix3_<Scalar> R  = Rz<Scalar>(x.get<0>().data());
  core::Vector3_<Scalar> x_ = x.get<1>().data()(0) * core::Vector3_<Scalar>(1.0, 0.0, 0.0);
  Scalar x_1_next           = (R * x_)(0);
  x_next.get<0>().setData(x.get<0>().data());
  x_next.get<1>().setData(core::Vector1_<Scalar>(x_1_next));
  return x_next;
}

bool string_in_array(const std::string& query, int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i)
    if (std::string(argv[i]) == query)
      return true;

  return false;
}

int main(int argc, char* argv[]) {
  JointType s0;
  CovType cov = CovType::Identity();

  s0.get<0>().setData(0.0);
  s0.get<1>().setData(core::Vector1_<Scalar>(1.0));

  SigmaPoints<JointType> spoints;

  JointType rec;
  CovType cov_rec;

  UnscentedTransform ut;

  if (string_in_array("HB", argc, argv)) {
    ut.weight_scheme_ = WeightScheme::HB;
  }

  ut.alpha_           = 4e-2;
  ut.cov_regularizer_ = 0;
  ut.toUnscented(s0, cov, spoints);

  for (int i = 0; i < spoints.points.size(); ++i) {
    std::cout << "-" << i << "-\n";
    const JointType& sp = spoints.points.at(i);
    std::cout << sp.get<0>().data() << "\n" << sp.get<1>().data().transpose() << "\n";
  }
  std::cout << "\n-applying f-\n\n";
  for (int i = 0; i < spoints.points.size(); ++i) {
    spoints.points.at(i) = f(spoints.points.at(i));
  }

  for (int i = 0; i < spoints.points.size(); ++i) {
    std::cout << "-" << i << "-\n";
    const JointType& sp = spoints.points.at(i);
    std::cout << sp.get<0>().data() << "\n" << sp.get<1>().data().transpose() << "\n";
  }

  ut.toMeanCov(spoints, rec, cov_rec);

  std::cout << "s0\n";
  std::cout << s0.get<0>().data() << "\n" << s0.get<1>().data() << "\n";
  std::cout << cov << "\n";
  std::cout << "rec\n";
  std::cout << rec.get<0>().data() << "\n" << rec.get<1>().data() << "\n";
  std::cout << cov_rec << std::endl;

  std::cout << "wm0: " << spoints.wm0 << "\n";
  std::cout << "wc0: " << spoints.wc0 << "\n";
  std::cout << "wmi: " << spoints.wmi << "\n";
  std::cout << "wci: " << spoints.wci << "\n";

  std::cout << "sum mean weights: " << spoints.wm0 + (spoints.points.size() - 1) * spoints.wmi
            << "\n";
  std::cout << "sum cov weights: " << spoints.wc0 + (spoints.points.size() - 1) * spoints.wci
            << "\n";
}