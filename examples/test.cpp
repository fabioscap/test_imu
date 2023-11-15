#include "common/common.h"
#include "common/manifold.h"
#include "common/manifold_impl.cpp"
#include "common/unscented.h"
#include "common/unscented_impl.cpp"
#include <random>
#include <srrg_geometry/geometry3d.h>

using namespace test_imu;

core::Vector3_<Scalar>
f(const core::Vector3_<Scalar>& x, const core::Matrix3_<Scalar>& A, const Scalar sigma) {
  // Create a random number generator with a normal distribution
  std::default_random_engine generator;
  std::normal_distribution<Scalar> distribution(0.0, sigma);

  // Generate random noise for each component
  Scalar noiseX = distribution(generator);
  Scalar noiseY = distribution(generator);
  Scalar noiseZ = distribution(generator);

  return A * x + core::Vector3_<Scalar>(noiseX, noiseY, noiseZ);
}

int main() {
  // linear transition function for euclidean state
  core::Matrix3_<Scalar> A = core::Matrix3_<Scalar>::Identity();
  A(1, 2)                  = 1.035;
  A(2, 2)                  = 0.3;
  A(2, 1)                  = 0.35;
  A(1, 2)                  = -0.35;

  float sigma = 0.;

  // test noise propagation with first order propagation
  size_t n_steps             = 10;
  core::Vector3_<Scalar> x_0 = core::Vector3_<Scalar>(1, -1, 0.5);

  core::Matrix3_<Scalar> cov_fo = core::Matrix3_<Scalar>::Identity();
  core::Vector3_<Scalar> x_fo   = x_0;

  for (size_t i = 0; i < n_steps; ++i) {
    cov_fo = A * cov_fo * A.transpose();
    x_fo   = f(x_fo, A, sigma);
  }

  std::cout << "FO\n";
  std::cout << "x: " << x_fo.transpose() << "\n";
  std::cout << "cov: \n" << cov_fo << "\n\n";

  core::Matrix3_<Scalar> cov_ut = core::Matrix3_<Scalar>::Identity();

  using StateType = Euclidean_<3>;

  StateType x_ut;
  x_ut.setData(x_0);

  SigmaPoints<StateType> spoints;

  UnscentedTransform ut;
  for (size_t i = 0; i < n_steps; ++i) {
    // get sigma points
    ut.toUnscented(x_ut, cov_ut, spoints);
    for (size_t j = 0; j < spoints.points.size(); ++j) {
      const auto& value = spoints.points.at(j).data();
      spoints.points.at(j).setData(f(value, A, sigma));
    }
    ut.toMeanCov(spoints, x_ut, cov_ut);
  }
  std::cout << "UT\n";
  std::cout << "x: " << x_fo.transpose() << "\n";
  std::cout << "cov: \n" << cov_fo << "\n\n";
}