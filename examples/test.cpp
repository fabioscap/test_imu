#include <srrg_geometry/geometry3d.h>

#include "synthetic/synthetic.h"

#include "variables_and_factors/imu_preintegration_factor.h"

#include <fstream>

int main() {
  using namespace test_imu;
  using TrajectoryType = SE3StraightTrajectory;
  float T              = 10;
  float freq           = 10;

  std::shared_ptr<TrajectoryType> traj = std::make_shared<TrajectoryType>(T);
  FakeImu imu(traj, freq, 102030);

  std::vector<std::pair<ImuMeasurement, srrg2_core::Isometry3f>> data;

  core::Vector3f pos, vel, acc;

  for (float t = 0; t < T; t += 0.1) {
    traj->sampleTrajectory(t, pos, vel, acc);
    std::cout << "pos: " << pos.transpose() << "\n";
    std::cout << "vel: " << vel.transpose() << "\n";
    std::cout << "acc: " << acc.transpose() << "\n";
    std::cout << "-\n";
  }
}