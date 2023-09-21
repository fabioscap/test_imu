#include <srrg_geometry/geometry3d.h>

#include "synthetic/synthetic.h"

#include "variables_and_factors/imu_preintegration_factor.h"

#include <fstream>

int main() {
  using namespace test_imu;
  using TrajectoryType = SE3CircleTrajectory;

  std::ofstream out_pred("/workspace/src/test_imu/examples/output_pred.txt");
  std::ofstream out_gt("/workspace/src/test_imu/examples/output_gt.txt");
  float T    = 10;
  float freq = 10;

  std::shared_ptr<TrajectoryType> traj = std::make_shared<TrajectoryType>(T);
  FakeImu imu(traj, freq, 102030);

  std::vector<std::pair<ImuMeasurement, srrg2_core::Isometry3f>> data;

  imu.generateData(data);

  srrg2_core::Isometry3f pose;
  srrg2_core::Vector3f vel;
  float t              = 0;
  ImuMeasurement& meas = data.at(0).first;

  // standard euler integration
  /*
  for (size_t i = 0; i < data.size(); ++i) {
    if (i == 0) {
      pose = data.at(i).second;
      t    = meas.timestamp;
      srrg2_core::Vector3f pos, acc;
      imu.trajectory().sampleTrajectory(t, pos, vel, acc);
      out_gt << data.at(i).second.matrix() << "\n";
      out_pred << pose.matrix() << "\n";
      continue;
    }
    out_gt << data.at(i).second.matrix() << "\n";
    out_pred << pose.matrix() << "\n";

    // srrg2_core::Isometry3f& T = data.at(i).second;
    ImuMeasurement& new_meas    = data.at(i).first;
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
  ImuPreintegrator integrator;

  core::Isometry3f initial_pose;
  for (size_t i = 0; i < data.size(); ++i) {
    if (i == 0) {
      initial_pose = data.at(i).second;
      meas         = data.at(i).first;
      integrator.reset(meas);
      srrg2_core::Vector3f pos, acc;
      imu.trajectory().sampleTrajectory(meas.timestamp, pos, vel, acc);
      out_gt << initial_pose.matrix() << "\n";
      out_pred << initial_pose.matrix() << "\n";
      continue;
    }

    ImuMeasurement new_meas = data.at(i).first;
    integrator.preintegrate(new_meas);

    t = new_meas.timestamp;
    srrg2_core::Vector3f posg, velg, accg;

    imu.trajectory().sampleTrajectory(meas.timestamp, posg, velg, accg);
    core::Isometry3f pose;
    core::Vector3f vel_now;

    integrator.getPrediction(initial_pose, vel, pose, vel_now);

    out_gt << data.at(i).second.matrix() << "\n";
    out_pred << pose.matrix() << "\n";

    meas = new_meas;
  }
}