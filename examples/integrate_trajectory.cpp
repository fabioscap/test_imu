#include <srrg_geometry/geometry3d.h>

#include "synthetic/synthetic.h"

#include "imu_preintegrator/imu_preintegrator.h"
#include "imu_preintegrator/imu_preintegrator_slim.h"
#include "imu_preintegrator/imu_preintegrator_ukf.h"

#include "common/manifold_impl.cpp"

#include <fstream>

int main() {
  using namespace test_imu;

  using PreintegratorType = ImuPreintegrator;
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

  // standard euler integration
  /*
  for (size_t i = 0; i < data.size(); ++i) {
    if (i == 0) {
      pose = std::get<2>(data.at(i));
      t    = meas.timestamp;
      srrg2_core::Vector3f pos, acc;
      imu.trajectory().sampleTrajectory(t, pos, vel, acc);
      out_gt << std::get<2>(data.at(i)).matrix() << "\n";
      out_pred << pose.matrix() << "\n";
      continue;
    }
    out_gt << std::get<2>(data.at(i)).matrix() << "\n";
    out_pred << pose.matrix() << "\n";

    // srrg2_core::Isometry3f& T = std::get<2>(data.at(i));
    ImuMeasurement& new_meas    = std::get<0>(data.at(i));
    float dt                    = new_meas.timestamp - t;
    Eigen::Vector3f translation = pose.translation();
    srrg2_core::Matrix3f R      = pose.rotation();

    srrg2_core::Vector3f posg, velg, accg;
    imu.trajectory().sampleTrajectory(t, posg, velg, accg);
    translation += vel * dt + 0.5 * R * meas.acceleration * dt * dt;
    vel += R * meas.acceleration * dt;
    R *=
      srrg2_core::geometry3d::expMapSO3(static_cast<srrg2_core::Vector3f>(dt * meas.angular_vel));
    srrg2_core::fixRotation(R);

    pose.translation() = translation;
    pose.linear()      = R;
    meas               = new_meas;
    t                  = new_meas.timestamp;
    std::cout << vel.x() << " " << vel.y() << "\n";
    std::cout << velg.x() << " " << velg.y() << "\n";
    std::cout << "---\n";
  }*/

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