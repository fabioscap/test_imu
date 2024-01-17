/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

// GTSAM related includes.
#include <gtsam/inference/Symbol.h>
#include <gtsam/navigation/CombinedImuFactor.h>
#include <gtsam/navigation/GPSFactor.h>
#include <gtsam/navigation/ImuFactor.h>
#include <gtsam/nonlinear/ISAM2.h>
#include <gtsam/nonlinear/ISAM2Params.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/slam/PriorFactor.h>
#include <gtsam/slam/dataset.h>

#include <gtsam/nonlinear/GaussNewtonOptimizer.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>

#include <cstring>
#include <fstream>
#include <iostream>

#include "common/common.h"

#include "synthetic/synthetic.h"

#include <chrono>

using namespace std;
using namespace gtsam;

using symbol_shorthand::B; // Bias  (ax,ay,az,gx,gy,gz)
using symbol_shorthand::V; // Vel   (xdot,ydot,zdot)
using symbol_shorthand::X; // Pose3 (x,y,z,r,p,y)

struct ImuMeasurement {
  double time;
  double dt;
  Vector3 accelerometer;
  Vector3 gyroscope; // omega
};

struct GpsMeasurement {
  double time;
  Vector3 position; // x,y,z
};

const string output_filename = "/workspace/src/test_imu/synthetic_gtsam.csv";

struct Sigmas {
  float acc       = 0.00175f;
  float gyro      = 0.00175f;
  float bias_acc  = 0.000167f;
  float bias_gyro = 0.000167f;

  float gps = 0.0002f;
} sigmas;

void generateData(std::vector<GpsMeasurement>& gps_measurements,
                  std::vector<ImuMeasurement>& imu_measurements,
                  Vector3& vel_zero) {
  using TrajectoryType = test_imu::SE3EightTrajectory;

  // total time of trajectory
  float T = 300;

  // imu measurements /sec
  float imu_freq = 50;

  // insert a gps measurement every gps_freq imu_measurements
  int gps_freq = 100;

  std::shared_ptr<TrajectoryType> traj = std::make_shared<TrajectoryType>(T);
  test_imu::FakeImu imu(traj, imu_freq, 102030);

  imu.std_acc()       = sigmas.acc;
  imu.std_gyro()      = sigmas.gyro;
  imu.std_bias_acc()  = sigmas.bias_acc;
  imu.std_bias_gyro() = sigmas.bias_gyro;

  std::vector<std::tuple<test_imu::ImuMeasurement, Eigen::Vector3f, Eigen::Isometry3f>> data;
  imu.generateData(data, true);

  std::random_device rd;
  std::mt19937 gen(rd());

  std::normal_distribution<double> distr(0, sigmas.gps);

  imu_measurements.clear();
  gps_measurements.clear();

  for (size_t i = 0; i < data.size(); ++i) {
    // add the imu measurements to the vector
    test_imu::ImuMeasurement meas_tmp = std::get<0>(data.at(i));
    ImuMeasurement meas;
    meas.dt = 1 / imu_freq;
    meas.accelerometer << meas_tmp.acceleration.cast<double>();
    meas.gyroscope << meas_tmp.angular_vel.cast<double>();
    meas.time = meas_tmp.timestamp;
    imu_measurements.push_back(meas);

    if (i % gps_freq == 0) {
      // add a new gps measurement

      Vector3 tr;
      tr << std::get<2>(data.at(i)).translation().cast<double>();
      GpsMeasurement gps_meas;
      // get the time from IMU
      gps_meas.time = std::get<0>(data.at(i)).timestamp;
      // get translation and add noise
      gps_meas.position = tr + Vector3(distr(gen), distr(gen), distr(gen));
      gps_measurements.push_back(gps_meas);
    }

    if (i == 0)
      vel_zero = std::get<1>(data.at(i)).cast<double>();
  }
}
int main(int argc, char* argv[]) {
  Values result;

  vector<ImuMeasurement> imu_measurements;
  vector<GpsMeasurement> gps_measurements;

  Vector3 vel0;

  generateData(gps_measurements, imu_measurements, vel0);
  std::cout << " loaded measurements\n";
  // Configure different variables
  // double t_offset = gps_measurements[0].time;
  size_t first_gps = 0;
  size_t gps_skip  = 1; // Skip this many GPS measurements each time
  double g         = 0.0;
  auto w_coriolis  = Vector3::Zero(); // zero vector

  LevenbergMarquardtParams parameters;
  parameters.setDiagonalDamping(false);
  parameters.setlambdaInitial(100);
  parameters.setlambdaFactor(1.0);
  parameters.setMaxIterations(50);

  parameters.setVerbosity("Summary");

  // Configure noise models
  auto noise_model_gps =
    noiseModel::Diagonal::Precisions((Vector6() << Vector3::Constant(0), sigmas.gps).finished());

  // Set initial conditions for the estimated trajectory
  // initial pose is the reference frame (navigation frame)
  auto current_pose_global = Pose3(Rot3(), gps_measurements[first_gps].position);
  // the vehicle is stationary at the beginning at position 0,0,0
  Vector3 current_velocity_global = vel0;
  auto current_bias               = imuBias::ConstantBias(); // init with zero bias

  auto sigma_init_x = noiseModel::Diagonal::Precisions(
    (Vector6() << Vector3::Constant(0), Vector3::Constant(1.0)).finished());
  auto sigma_init_v = noiseModel::Diagonal::Sigmas(Vector3::Constant(1000.0));
  auto sigma_init_b = noiseModel::Diagonal::Sigmas(
    (Vector6() << Vector3::Constant(0.100), Vector3::Constant(5.00e-05)).finished());

  // Set IMU preintegration parameters
  Matrix33 measured_acc_cov   = I_3x3 * sigmas.acc;
  Matrix33 measured_omega_cov = I_3x3 * sigmas.gyro;
  // error committed in integrating position from velocities
  Matrix33 integration_error_cov = I_3x3 * 0.0;

  auto imu_params                     = PreintegratedImuMeasurements::Params::MakeSharedU(g);
  imu_params->accelerometerCovariance = measured_acc_cov;      // acc white noise in continuous
  imu_params->integrationCovariance   = integration_error_cov; // integration uncertainty continuous
  imu_params->gyroscopeCovariance     = measured_omega_cov;    // gyro white noise in continuous
  imu_params->omegaCoriolis           = w_coriolis;

  std::shared_ptr<PreintegratedImuMeasurements> current_summarized_measurement = nullptr;

  NonlinearFactorGraph graph;
  Values values;

  size_t j                              = 0;
  size_t included_imu_measurement_count = 0;

  for (size_t i = first_gps; i < gps_measurements.size(); i = i + gps_skip) {
    std::cout << std::fixed << "gps_idx: " << i << " t: " << gps_measurements[i].time << "\n";

    // At each non=IMU measurement we initialize a new node in the graph
    auto current_pose_key = X(i);
    auto current_vel_key  = V(i);
    auto current_bias_key = B(i);

    if (i == first_gps) {
      // Create initial estimate and prior on initial pose, velocity, and biases
      values.insert(current_pose_key, current_pose_global);
      values.insert(current_vel_key, current_velocity_global);
      values.insert(current_bias_key, current_bias);
      graph.emplace_shared<PriorFactor<Pose3>>(current_pose_key, current_pose_global, sigma_init_x);
      graph.emplace_shared<PriorFactor<Vector3>>(
        current_vel_key, current_velocity_global, sigma_init_v);
      graph.emplace_shared<PriorFactor<imuBias::ConstantBias>>(
        current_bias_key, current_bias, sigma_init_b);
      continue;
    }
    int prev_idx            = (i == first_gps + 1) ? first_gps : i - gps_skip;
    const double t_previous = gps_measurements.at(prev_idx).time;
    double t                = gps_measurements[i].time;

    // Summarize IMU data between the previous GPS measurement and now
    current_summarized_measurement =
      std::make_shared<PreintegratedImuMeasurements>(imu_params, current_bias);

    auto start = std::chrono::high_resolution_clock::now();
    while (j < imu_measurements.size() && imu_measurements[j].time <= t) {
      if (imu_measurements[j].time >= t_previous) {
        current_summarized_measurement->integrateMeasurement(imu_measurements[j].accelerometer,
                                                             imu_measurements[j].gyroscope,
                                                             imu_measurements.at(j + 1).time -
                                                               imu_measurements.at(j).time);
      }
      j++;
    }
    auto stop     = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

    std::cout << "Time taken: " << duration.count() << " nanoseconds" << std::endl;

    // Create IMU factor
    auto previous_pose_key = X(prev_idx);
    auto previous_vel_key  = V(prev_idx);
    auto previous_bias_key = B(prev_idx);

    graph.emplace_shared<ImuFactor>(previous_pose_key,
                                    previous_vel_key,
                                    current_pose_key,
                                    current_vel_key,
                                    previous_bias_key,
                                    *current_summarized_measurement);

    // Bias evolution as given in the IMU metadata
    auto sigma_between_b = noiseModel::Diagonal::Sigmas(
      (Vector6() << Vector3::Constant(sqrt(included_imu_measurement_count) *
                                      std::sqrt(sigmas.bias_acc)),
       Vector3::Constant(sqrt(included_imu_measurement_count) * std::sqrt(sigmas.bias_gyro)))
        .finished());
    graph.emplace_shared<BetweenFactor<imuBias::ConstantBias>>(
      previous_bias_key, current_bias_key, imuBias::ConstantBias(), sigma_between_b);

    // Create GPS factor
    auto gps_pose = Pose3(current_pose_global.rotation(), gps_measurements[i].position);

    graph.emplace_shared<PriorFactor<Pose3>>(current_pose_key, gps_pose, noise_model_gps);

    printf("############ POSE INCLUDED AT TIME %.6lf ############\n", t);
    cout << gps_pose.translation();
    printf("\n\n");

    // Add initial values for velocity and bias based on the previous
    // estimates
    values.insert(current_pose_key, gps_pose);
    values.insert(current_vel_key, current_velocity_global);
    values.insert(current_bias_key, current_bias);

    // Stop iterating once the change in error between steps is less than this value
    parameters.relativeErrorTol = 1e-5;
    // Do not perform more than N ite

    LevenbergMarquardtOptimizer optimizer(graph, values, parameters);
    Values optimizedValues = optimizer.optimize();
    values                 = optimizedValues;

    current_pose_global     = values.at<Pose3>(current_pose_key);
    current_velocity_global = values.at<Vector3>(current_vel_key);
    current_bias            = values.at<imuBias::ConstantBias>(current_bias_key);
  }

  // Save results to file
  printf("\nWriting results to file...\n");
  FILE* fp_out = fopen(output_filename.c_str(), "w+");

  graph.print();

  for (size_t i = 0; i < gps_measurements.size(); i = i + gps_skip) {
    auto pose_key = X(i);
    auto vel_key  = V(i);
    auto bias_key = B(i);

    auto pose     = values.at<Pose3>(pose_key);
    auto velocity = values.at<Vector3>(vel_key);
    auto bias     = values.at<imuBias::ConstantBias>(bias_key);

    auto pose_euler = pose.rotation().rpy();

    using namespace test_imu;

    std::cout << pose.rotation() << "\n"
              << Rx(pose_euler.x()) * Ry(pose_euler.y()) * Rz(pose_euler.z()) << "\n\n";

    cout << "State at #" << i << endl;
    cout << "Pose:" << endl << pose << endl;
    cout << "Velocity:" << endl << velocity << endl;
    cout << "Bias:" << endl << bias << endl;

    fprintf(fp_out,
            "%f,%f,%f,%f,%f,%f,%f\n",
            gps_measurements[i].time,
            pose.x(),
            pose.y(),
            pose.z(),
            pose_euler.x(),
            pose_euler.y(),
            pose_euler.z());
  }

  fclose(fp_out);
}