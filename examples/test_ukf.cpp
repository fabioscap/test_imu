#include "common/manifold.h"
#include "common/unscented.h"
#include <srrg_geometry/geometry3d.h>

using namespace test_imu;

int main() {
  using StateType = ManifoldSO3;
  using InputType = Euclidean_<3>;
  using JointType = ManifoldComp_<StateType, InputType>;

  using CovType = UnscentedTransform::CovType<JointType>;

  JointType s0;
}