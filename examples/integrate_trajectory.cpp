#include <srrg_geometry/geometry3d.h>

#include "synthetic/synthetic.h"

#include "imu_preintegrator/imu_preintegrator.h"
#include "imu_preintegrator/imu_preintegrator_ukf.h"

#include "common/manifold_impl.cpp"

#include <fstream>

int main() {
  using namespace test_imu;

  using PreintegratorType = ImuPreintegratorUKFSlim;
  using TrajectoryType    = SE3EightTrajectory;

  std::ofstream out_pred("/workspace/src/test_imu/examples/output_pred.txt");
  std::ofstream out_gt("/workspace/src/test_imu/examples/output_gt.txt");
  float T    = 10;
  float freq = 50;

  std::shared_ptr<TrajectoryType> traj = std::make_shared<TrajectoryType>(T);
  FakeImu imu(traj, freq, 102030);

  std::vector<std::tuple<ImuMeasurement, core::Vector3f, core::Isometry3f>> data;

  imu.generateData(data);

  srrg2_core::Isometry3f pose;
  srrg2_core::Vector3f vel;
  float t              = 0;
  ImuMeasurement& meas = std::get<0>(data.at(0));

  // imu preintegration
  PreintegratorType integrator;
  integrator.reset();
  float dt = 1 / imu.freq();
  std::cout << "dt: " << dt << std::endl;

  core::Isometry3f initial_pose;
  for (size_t i = 0; i < data.size(); ++i) {
    if (i == 0) {
      initial_pose = std::get<2>(data.at(i));
      vel          = std::get<1>(data.at(i));
      meas         = std::get<0>(data.at(i));
      integrator.preintegrate(meas, dt);
      srrg2_core::Vector3f pos, acc;

      out_gt << initial_pose.matrix() << "\n";
      out_pred << initial_pose.matrix() << "\n";
      continue;
    }

    ImuMeasurement new_meas = std::get<0>(data.at(i));
    integrator.preintegrate(new_meas, dt);

    t = new_meas.timestamp;

    core::Isometry3f pose;
    core::Vector3f vel_now;

    integrator.getPrediction(initial_pose, vel, pose, vel_now);

    out_gt << std::get<2>(data.at(i)).matrix() << "\n";
    out_pred << pose.matrix() << "\n";

    meas = new_meas;
  }
}