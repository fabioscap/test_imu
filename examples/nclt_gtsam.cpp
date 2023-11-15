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

#include <cstring>
#include <fstream>
#include <iostream>

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

const string output_filename = "/workspace/src/test_imu/nclt_gtsam.csv";

struct Sigmas {
  float acc       = 0.00001f;
  float gyro      = 0.00001f;
  float bias_acc  = 0.0000001f;
  float bias_gyro = 0.0000001f;

  float gps = 0.02f;
} sigmas;

void loadNCLTData(vector<ImuMeasurement>& imu_measurements,
                  vector<GpsMeasurement>& gps_measurements) {
  const std::string data_folder("/workspace/src/test_imu/data/nclt");
  string imu_filename = "ms25.csv";
  string gps_filename = "gps_rtk.csv";
  ifstream imu_file(data_folder + "/" + imu_filename);
  if (!imu_file) {
    throw std::runtime_error("loadNCLTData: cannot open IMU data file at" + data_folder + "/" +
                             imu_filename + ".");
  }

  string line;
  while (std::getline(imu_file, line)) {
    ImuMeasurement meas;

    sscanf(line.c_str(),
           "%lf,%*f,%*f,%*f,%lf,%lf,%lf,%lf,%lf,%lf",
           &meas.time,
           &meas.accelerometer.x(),
           &meas.accelerometer.y(),
           &meas.accelerometer.z(),
           &meas.gyroscope.x(),
           &meas.gyroscope.y(),
           &meas.gyroscope.z());
    meas.time *= 1e-6;
    if (imu_measurements.size() == 0) {
    } else if (imu_measurements.size() == 1) {
      meas.dt                   = meas.time - imu_measurements.at(0).time;
      imu_measurements.at(0).dt = meas.dt;
    } else {
      meas.dt = meas.time - imu_measurements.at(imu_measurements.size() - 2).time;
    }

    imu_measurements.push_back(meas);
  }

  ifstream gps_file(data_folder + "/" + gps_filename);
  if (!gps_file) {
    throw std::runtime_error("loadNCLTData: cannot open GPS data file at" + data_folder + "/" +
                             gps_filename + ".");
  }

  double lat0 = (42.293227f / 180.0f) * M_PI;
  double lng0 = (-83.709657 / 180.0f) * M_PI;
  double alt0 = 270.0;

  double re = 6378135;
  double rp = 6356750;

  double den = pow(re * cos(lat0), 2) + pow(rp * sin(lat0), 2.0);

  // radius north south
  double rns = pow(re * rp, 2) / pow(den, 3.0 / 2.0);

  // radius east west
  double rew = (re * re) / sqrt(den);

  while (std::getline(gps_file, line)) {
    GpsMeasurement meas;

    double lat, lng, alt;

    int fix_mode;

    sscanf(line.c_str(), "%lf,%d,%*d,%lf,%lf,%lf,%*f,%*f", &meas.time, &fix_mode, &lat, &lng, &alt);

    if (fix_mode == 3) {
      // convert using formulas in pdf
      meas.position.x() = sin(lat - lat0) * rns;
      meas.position.y() = sin(lng - lng0) * rew * cos(lat0);
      meas.position.z() = alt0 - alt;
      meas.time *= 1e-6;
      gps_measurements.push_back(meas);
    }
  }
}

int main(int argc, char* argv[]) {
  // DOES NOT WORK: ILL CONDITIONED SYSTEM
  Values result;

  vector<ImuMeasurement> imu_measurements;
  vector<GpsMeasurement> gps_measurements;
  loadNCLTData(imu_measurements, gps_measurements);
  std::cout << " loaded measurements\n";
  // Configure different variables
  // double t_offset = gps_measurements[0].time;
  size_t first_gps_pose = 1;
  size_t gps_skip       = 10; // Skip this many GPS measurements each time
  double g              = 9.8;
  auto w_coriolis       = Vector3::Zero(); // zero vector

  // Configure noise models
  auto noise_model_gps = noiseModel::Diagonal::Precisions(
    (Vector6() << Vector3::Constant(0), Vector3::Constant(1.0 / sigmas.gps)).finished());

  // Set initial conditions for the estimated trajectory
  // initial pose is the reference frame (navigation frame)
  auto current_pose_global = Pose3(Rot3(), gps_measurements[first_gps_pose].position);
  // the vehicle is stationary at the beginning at position 0,0,0
  Vector3 current_velocity_global = Vector3::Zero();
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

  // Set ISAM2 parameters and create ISAM2 solver object
  ISAM2Params isam_params;
  isam_params.factorization   = ISAM2Params::CHOLESKY;
  isam_params.relinearizeSkip = 10;

  ISAM2 isam(isam_params);

  // Create the factor graph and values object that will store new factors and
  // values to add to the incremental graph
  NonlinearFactorGraph new_factors;
  Values new_values; // values storing the initial estimates of new nodes in
                     // the factor graph

  printf("-- Starting main loop: inference is performed at each time step, but we "
         "plot trajectory every 10 steps\n");
  size_t j                              = 0;
  size_t included_imu_measurement_count = 0;

  for (size_t i = first_gps_pose; i < gps_measurements.size() - 1; i++) {
    // At each non=IMU measurement we initialize a new node in the graph
    auto current_pose_key = X(i);
    auto current_vel_key  = V(i);
    auto current_bias_key = B(i);
    double t              = gps_measurements[i].time;

    if (i == first_gps_pose) {
      // Create initial estimate and prior on initial pose, velocity, and biases
      new_values.insert(current_pose_key, current_pose_global);
      new_values.insert(current_vel_key, current_velocity_global);
      new_values.insert(current_bias_key, current_bias);
      new_factors.emplace_shared<PriorFactor<Pose3>>(
        current_pose_key, current_pose_global, sigma_init_x);
      new_factors.emplace_shared<PriorFactor<Vector3>>(
        current_vel_key, current_velocity_global, sigma_init_v);
      new_factors.emplace_shared<PriorFactor<imuBias::ConstantBias>>(
        current_bias_key, current_bias, sigma_init_b);
    } else {
      double t_previous = gps_measurements[i - 1].time;

      // Summarize IMU data between the previous GPS measurement and now
      current_summarized_measurement =
        std::make_shared<PreintegratedImuMeasurements>(imu_params, current_bias);

      while (j < imu_measurements.size() && imu_measurements[j].time <= t) {
        if (imu_measurements[j].time >= t_previous) {
          current_summarized_measurement->integrateMeasurement(imu_measurements[j].accelerometer,
                                                               imu_measurements[j].gyroscope,
                                                               imu_measurements[j].dt);
          included_imu_measurement_count++;
        }
        j++;
      }

      // Create IMU factor
      auto previous_pose_key = X(i - 1);
      auto previous_vel_key  = V(i - 1);
      auto previous_bias_key = B(i - 1);

      new_factors.emplace_shared<ImuFactor>(previous_pose_key,
                                            previous_vel_key,
                                            current_pose_key,
                                            current_vel_key,
                                            previous_bias_key,
                                            *current_summarized_measurement);
      std::cout << "cov:\n" << current_summarized_measurement->preintMeasCov() << "\n\n";
      std::cout << "determinant:\n"
                << current_summarized_measurement->preintMeasCov().determinant() << "\n\n";
      std::cout << "omega:\n"
                << current_summarized_measurement->preintMeasCov().inverse() << "\n\n";

      // Bias evolution as given in the IMU metadata
      auto sigma_between_b = noiseModel::Diagonal::Sigmas(
        (Vector6() << Vector3::Constant(sqrt(included_imu_measurement_count) *
                                        std::sqrt(sigmas.bias_acc)),
         Vector3::Constant(sqrt(included_imu_measurement_count) * std::sqrt(sigmas.bias_gyro)))
          .finished());
      new_factors.emplace_shared<BetweenFactor<imuBias::ConstantBias>>(
        previous_bias_key, current_bias_key, imuBias::ConstantBias(), sigma_between_b);

      // Create GPS factor
      auto gps_pose = Pose3(current_pose_global.rotation(), gps_measurements[i].position);
      if ((i % gps_skip) == 0) {
        new_factors.emplace_shared<PriorFactor<Pose3>>(current_pose_key, gps_pose, noise_model_gps);
        new_values.insert(current_pose_key, gps_pose);

        printf("############ POSE INCLUDED AT TIME %.6lf ############\n", t);
        cout << gps_pose.translation();
        printf("\n\n");
      } else {
        new_values.insert(current_pose_key, current_pose_global);
      }

      // Add initial values for velocity and bias based on the previous
      // estimates
      new_values.insert(current_vel_key, current_velocity_global);
      new_values.insert(current_bias_key, current_bias);

      // Update solver
      // =======================================================================
      // We accumulate 2*GPSskip GPS measurements before updating the solver at
      // first so that the heading becomes observable.
      if (i > (first_gps_pose + 2 * gps_skip)) {
        printf("############ NEW FACTORS AT TIME %.6lf ############\n", t);
        new_factors.print();

        isam.update(new_factors, new_values);

        // Reset the newFactors and newValues list
        new_factors.resize(0);
        new_values.clear();

        // Extract the result/current estimates
        result = isam.calculateEstimate();

        current_pose_global     = result.at<Pose3>(current_pose_key);
        current_velocity_global = result.at<Vector3>(current_vel_key);
        current_bias            = result.at<imuBias::ConstantBias>(current_bias_key);

        printf("\n############ POSE AT TIME %lf ############\n", t);
        current_pose_global.print();
        printf("\n\n");
      }
    }
  }

  // Save results to file
  printf("\nWriting results to file...\n");
  FILE* fp_out = fopen(output_filename.c_str(), "w+");
  fprintf(fp_out, "#time(s),x(m),y(m),z(m),qx,qy,qz,qw,gt_x(m),gt_y(m),gt_z(m)\n");

  result = isam.calculateEstimate();
  for (size_t i = first_gps_pose; i < gps_measurements.size() - 1; i++) {
    auto pose_key = X(i);
    auto vel_key  = V(i);
    auto bias_key = B(i);

    auto pose     = result.at<Pose3>(pose_key);
    auto velocity = result.at<Vector3>(vel_key);
    auto bias     = result.at<imuBias::ConstantBias>(bias_key);

    auto pose_quat = pose.rotation().toQuaternion();
    auto gps       = gps_measurements[i].position;

    cout << "State at #" << i << endl;
    cout << "Pose:" << endl << pose << endl;
    cout << "Velocity:" << endl << velocity << endl;
    cout << "Bias:" << endl << bias << endl;

    fprintf(fp_out,
            "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
            gps_measurements[i].time,
            pose.x(),
            pose.y(),
            pose.z(),
            pose_quat.x(),
            pose_quat.y(),
            pose_quat.z(),
            pose_quat.w(),
            gps(0),
            gps(1),
            gps(2));
  }

  fclose(fp_out);
}