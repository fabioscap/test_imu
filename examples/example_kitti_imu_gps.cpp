#include <fstream>
#include <iostream>

#define _USE_MATH_DEFINES
#include "variables_and_factors/gps_factor_ad.h"
#include "variables_and_factors/imu_preintegration_factor.h"

#include "variables_and_factors/instances.h"

#include "imu_preintegrator/imu_preintegrator.h"
#include "imu_preintegrator/imu_preintegrator_ukf.h"

#include "srrg_solver/solver_core/factor_graph.h"
#include "srrg_solver/solver_core/instances.h"

#include "srrg_solver/solver_core/internals/linear_solvers/instances.h"
#include "srrg_solver/solver_core/solver.h"
#include "srrg_solver/variables_and_factors/types_3d/all_types.h"
#include "srrg_solver/variables_and_factors/types_3d/instances.h"
#include <srrg_solver/solver_core/iteration_algorithm_ddl.h>
#include <srrg_solver/solver_core/iteration_algorithm_dl.h>
#include <srrg_solver/solver_core/iteration_algorithm_gn.h>
#include <srrg_solver/solver_core/iteration_algorithm_lm.h>
#include <srrg_system_utils/parse_command_line.h>

#include "srrg_solver/solver_core/factor_graph_view.h"

#include <thread>
#include <unistd.h>

#include <string.h>

using namespace std;
using namespace srrg2_core;
using namespace srrg2_solver;

const std::string example_folder("/workspace/src/test_imu");

struct KittiCalibration {
  double body_ptx;
  double body_pty;
  double body_ptz;
  double body_prx;
  double body_pry;
  double body_prz;
  double accelerometer_sigma;
  double gyroscope_sigma;
  double integration_sigma;
  double accelerometer_bias_sigma;
  double gyroscope_bias_sigma;
  double average_delta_t;
};

struct ImuMeasurement {
  double time;
  double dt;
  Vector3f accelerometer;
  Vector3f gyroscope; // omega
};

struct GpsMeasurement {
  double time;
  Vector3f position; // x,y,z
};

void loadKittiData(KittiCalibration& kitti_calibration,
                   vector<ImuMeasurement>& imu_measurements,
                   vector<GpsMeasurement>& gps_measurements);

void viewGraph(ViewerCanvasPtr canvas_);

bool string_in_array(const string& query, int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i)
    if (std::string(argv[i]) == query)
      return true;

  return false;
}

int main(int argc, char* argv[]) {
  // incremental optimization parameters
  const int window_size = 8;
  const int gps_skip    = 1; // 1 is no skips
  // which variables get locally optimized
  IdSet variable_ids;
  std::queue<VariableBase::Id> pose_local_ids;
  std::queue<VariableBase::Id> vel_local_ids;
  std::queue<VariableBase::Id> acc_bias_local_ids;
  std::queue<VariableBase::Id> gyro_bias_local_ids;

  bool use_ukf  = string_in_array("ukf", argc, argv);
  bool use_imu  = !string_in_array("noimu", argc, argv);
  bool use_slim = string_in_array("slim", argc, argv);
  bool use_gps  = !string_in_array("nogps", argc, argv);

  // inspect covariances
  ofstream cov_dump(example_folder + "/covariance.txt");
  ofstream det_dump(example_folder + "/determinant.txt");

  variables_and_factors_3d_registerTypes();
  variables_and_factors_imu_registerTypes();

  KittiCalibration kitti_calibration;
  vector<ImuMeasurement> imu_measurements;
  vector<GpsMeasurement> gps_measurements;
  loadKittiData(kitti_calibration, imu_measurements, gps_measurements);

  Solver solver;
  solver.param_termination_criteria.setValue(nullptr);
  solver.param_max_iterations.pushBack(50);
  IterationAlgorithmBasePtr alg(new IterationAlgorithmGN);
  std::shared_ptr<IterationAlgorithmGN> temp = std::dynamic_pointer_cast<IterationAlgorithmGN>(alg);
  temp->param_damping.setValue(1e3);
  solver.param_algorithm.setValue(alg);
  solver.param_verbose.setValue(true);

  FactorGraphPtr graph(new FactorGraph);

  using VarPoseImuType = VariableSE3ExpMapRightAD;
  using VarVelImuType  = VariableVector3AD;
  using ImuBiasVar     = VariableVector3AD;
  using FactorGpsType  = GpsFactorAD;
  using FactorBiasType = BiasErrorFactorAD;

  // initialization

  const Vector3f& init_gps_pose = gps_measurements.at(0).position;
  Isometry3f init_pose          = Isometry3f::Identity();
  init_pose.translation()       = init_gps_pose;
  init_pose.linear()            = test_imu::Rz(0.785f);

  VarPoseImuType* prev_pose_var = new VarPoseImuType();
  prev_pose_var->setEstimate(init_pose);
  prev_pose_var->setGraphId(0);
  graph->addVariable(VariableBasePtr(prev_pose_var));

  FactorGpsType* gps_factor = new FactorGpsType();
  gps_factor->setVariableId(0, 0);

  float info_gps = 1 / 0.007;
  gps_factor->setInformationMatrix(Matrix3f::Identity() * info_gps);
  gps_factor->setMeasurement(init_gps_pose);
  graph->addFactor(FactorBasePtr(gps_factor));

  // std::cout << "adding pose variable with id: " << prev_pose_var->graphId() << "\n";
  pose_local_ids.push(prev_pose_var->graphId());
  variable_ids.insert(prev_pose_var->graphId());

  VarVelImuType* prev_vel_var = new VarVelImuType();
  prev_vel_var->setEstimate(Vector3f::Zero());
  prev_vel_var->setGraphId(1);
  // prev_vel_var->setStatus(VariableBase::Status::Fixed);
  graph->addVariable(VariableBasePtr(prev_vel_var));
  vel_local_ids.push(prev_vel_var->graphId());
  variable_ids.insert(prev_vel_var->graphId());
  // std::cout << "adding vel variable with id: " << prev_vel_var->graphId() << "\n";

  ImuBiasVar* prev_bias_acc = new ImuBiasVar();
  prev_bias_acc->setGraphId(2);
  prev_bias_acc->setEstimate(Vector3f::Zero());
  prev_bias_acc->setStatus(VariableBase::Fixed);
  graph->addVariable(VariableBasePtr(prev_bias_acc));
  acc_bias_local_ids.push(prev_bias_acc->graphId());
  variable_ids.insert(prev_bias_acc->graphId());

  std::cout << "adding bias_acc variable with id: " << prev_bias_acc->graphId() << "\n";

  ImuBiasVar* prev_bias_gyro = new ImuBiasVar();
  prev_bias_gyro->setGraphId(3);
  prev_bias_gyro->setEstimate(Vector3f::Zero());
  prev_bias_gyro->setStatus(VariableBase::Fixed);
  graph->addVariable(VariableBasePtr(prev_bias_gyro));
  gyro_bias_local_ids.push(prev_bias_gyro->graphId());
  variable_ids.insert(prev_bias_gyro->graphId());

  std::cout << "adding bias_gyro variable with id: " << prev_bias_gyro->graphId() << "\n";

  size_t graph_id = 4;

  test_imu::ImuPreintegratorBase* imu_preintegrator;
  if (use_slim) {
    std::cout << "slim ";
    if (use_ukf) {
      std::cout << "ukf\n";
      imu_preintegrator = new test_imu::ImuPreintegratorUKFSlim();
    } else {
      std::cout << "normal\n";
      imu_preintegrator = new test_imu::ImuPreintegratorSlim();
    }
  } else {
    std::cout << "full ";
    if (use_ukf) {
      std::cout << "ukf\n";
      imu_preintegrator = new test_imu::ImuPreintegratorUKF();
    } else {
      std::cout << "normal\n";
      imu_preintegrator = new test_imu::ImuPreintegrator();
    }
  }

  imu_preintegrator->setNoiseGyroscope(Vector3f::Constant(kitti_calibration.gyroscope_sigma));
  imu_preintegrator->setNoiseAccelerometer(
    Vector3f::Constant(kitti_calibration.accelerometer_sigma));
  imu_preintegrator->setNoiseBiasGyroscope(
    Vector3f::Constant(kitti_calibration.gyroscope_bias_sigma));
  imu_preintegrator->setNoiseBiasAccelerometer(
    Vector3f::Constant(kitti_calibration.accelerometer_bias_sigma));

  imu_preintegrator->setNoiseGyroscope(
    Vector3f::Constant(kitti_calibration.gyroscope_sigma * kitti_calibration.gyroscope_sigma));
  imu_preintegrator->setNoiseAccelerometer(Vector3f::Constant(
    kitti_calibration.accelerometer_sigma * kitti_calibration.accelerometer_sigma));
  imu_preintegrator->setNoiseBiasGyroscope(Vector3f::Constant(
    kitti_calibration.gyroscope_bias_sigma * kitti_calibration.gyroscope_bias_sigma));
  imu_preintegrator->setNoiseBiasAccelerometer(Vector3f::Constant(
    kitti_calibration.accelerometer_bias_sigma * kitti_calibration.accelerometer_bias_sigma));

  size_t j = 0;
  for (size_t i = 1; i < gps_measurements.size(); i = i + gps_skip) {
    int prev_idx = (i == 1) ? 0 : i - gps_skip;

    const double t_previous = gps_measurements.at(prev_idx).time;
    const double gps_time   = gps_measurements.at(i).time;
    const double dT         = gps_time - t_previous;

    imu_preintegrator->reset();
    imu_preintegrator->setBiasAcc(prev_bias_acc->estimate());
    imu_preintegrator->setBiasGyro(prev_bias_gyro->estimate());

    while (j < imu_measurements.size() && imu_measurements.at(j).time <= gps_time) {
      if (imu_measurements.at(j).time >= t_previous) {
        test_imu::ImuMeasurement meas;
        meas.acceleration = imu_measurements.at(j).accelerometer;
        meas.angular_vel  = imu_measurements.at(j).gyroscope;
        meas.timestamp    = imu_measurements.at(j).time;
        imu_preintegrator->preintegrate(meas, imu_measurements.at(j).dt);
      }

      j++;
    }

    // order of var indeces
    // 1 pose from imu
    // 2 vel from imu
    // 3 pose to imu
    // 4 vel to imu

    const Vector3f curr_gps_pose = gps_measurements.at(i).position;
    Isometry3f curr_pose         = Isometry3f::Identity();
    curr_pose.translation()      = curr_gps_pose;
    curr_pose.linear()           = prev_pose_var->estimate().linear();

    VarPoseImuType* curr_pose_var = new VarPoseImuType();
    curr_pose_var->setEstimate(curr_pose);
    curr_pose_var->setGraphId(graph_id++);
    graph->addVariable(VariableBasePtr(curr_pose_var));
    pose_local_ids.push(curr_pose_var->graphId());
    variable_ids.insert(curr_pose_var->graphId());
    // std::cout << "adding pose variable with id: " << curr_pose_var->graphId() << "\n";

    VarVelImuType* curr_vel_var = new VarVelImuType();
    curr_vel_var->setGraphId(graph_id++);
    const Vector3f prev_vel =
      static_cast<VarVelImuType*>(graph->variable(graph_id - 3))->estimate();
    curr_vel_var->setEstimate(prev_vel);
    graph->addVariable(VariableBasePtr(curr_vel_var));
    vel_local_ids.push(curr_vel_var->graphId());
    variable_ids.insert(curr_vel_var->graphId());
    // std::cout << "adding vel variable with id: " << curr_vel_var->graphId() << "\n";

    ImuBiasVar* curr_bias_acc  = new ImuBiasVar();
    ImuBiasVar* curr_bias_gyro = new ImuBiasVar();
    curr_bias_acc->setGraphId(graph_id++);
    curr_bias_gyro->setGraphId(graph_id++);
    curr_bias_acc->setEstimate(Vector3f::Zero());
    curr_bias_gyro->setEstimate(Vector3f::Zero());
    graph->addVariable(VariableBasePtr(curr_bias_acc));
    acc_bias_local_ids.push(curr_bias_acc->graphId());
    variable_ids.insert(curr_bias_acc->graphId());

    graph->addVariable(VariableBasePtr(curr_bias_gyro));
    gyro_bias_local_ids.push(curr_bias_gyro->graphId());
    variable_ids.insert(curr_bias_gyro->graphId());

    std::cout << "imu preintegrator has absorbed: " << imu_preintegrator->measurements().size()
              << " measurements.\n";
    if (use_imu && imu_preintegrator->measurements().size() > 0) {
      if (!use_slim) {
        std::cout << "full factor\n";
        ImuPreintegrationFactorAD* imu_factor = new ImuPreintegrationFactorAD();
        imu_factor->grav(Vector3f(0.f, 0.f, -9.80655));
        imu_factor->setVariableId(0, prev_pose_var->graphId());
        imu_factor->setVariableId(1, prev_vel_var->graphId());
        imu_factor->setVariableId(2, curr_pose_var->graphId());
        imu_factor->setVariableId(3, curr_vel_var->graphId());

        imu_factor->setVariableId(4, prev_bias_acc->graphId());
        imu_factor->setVariableId(5, prev_bias_gyro->graphId());

        imu_factor->setVariableId(6, curr_bias_acc->graphId());
        imu_factor->setVariableId(7, curr_bias_gyro->graphId());
        imu_factor->setMeasurement(*imu_preintegrator);

        graph->addFactor(FactorBasePtr(imu_factor));
      } else {
        std::cout << "slim factor\n";
        ImuPreintegrationFactorSlimAD* imu_factor = new ImuPreintegrationFactorSlimAD();
        imu_factor->grav(Vector3f(0.f, 0.f, -9.80655));
        imu_factor->setVariableId(0, prev_pose_var->graphId());
        imu_factor->setVariableId(1, prev_vel_var->graphId());
        imu_factor->setVariableId(2, curr_pose_var->graphId());
        imu_factor->setVariableId(3, curr_vel_var->graphId());

        FactorBiasType* bias_factor = new FactorBiasType();
        bias_factor->setVariableId(0, prev_bias_acc->graphId());
        bias_factor->setVariableId(1, prev_bias_gyro->graphId());

        bias_factor->setVariableId(2, curr_bias_acc->graphId());
        bias_factor->setVariableId(3, curr_bias_gyro->graphId());
        Matrix6f bias_sigma = Matrix6f::Identity();
        bias_sigma.block<3, 3>(0, 0) *= dT * kitti_calibration.accelerometer_bias_sigma;
        bias_sigma.block<3, 3>(3, 3) *= dT * kitti_calibration.gyroscope_bias_sigma;
        bias_factor->setInformationMatrix(bias_sigma.inverse());
        imu_factor->setMeasurement(*imu_preintegrator);

        graph->addFactor(FactorBasePtr(bias_factor));
        graph->addFactor(FactorBasePtr(imu_factor));
      }
    }

    FactorGpsType* gps_factor = new FactorGpsType();
    gps_factor->setVariableId(0, curr_pose_var->graphId());

    gps_factor->setInformationMatrix(Matrix3f::Identity() * info_gps);
    gps_factor->setMeasurement(curr_gps_pose);

    if (use_gps) {
      graph->addFactor(FactorBasePtr(gps_factor));
    }
    cov_dump << "cov:\n" << imu_preintegrator->sigma() << "\n\n";
    cov_dump << "omega:\n" << imu_preintegrator->sigma().inverse() << "\n\n";
    det_dump << imu_preintegrator->sigma().determinant() << "\n";
    imu_preintegrator->reset();

    // local optimization
    std::cout << "#local poses " << pose_local_ids.size() << "\n";
    FactorGraphView local_window;

    std::cout << "local opt.\n";

    while (pose_local_ids.size() > window_size) {
      std::cout << "sliding variables\n";
      variable_ids.erase(pose_local_ids.front());
      variable_ids.erase(vel_local_ids.front());
      variable_ids.erase(acc_bias_local_ids.front());
      variable_ids.erase(gyro_bias_local_ids.front());
      pose_local_ids.pop();
      vel_local_ids.pop();
      acc_bias_local_ids.pop();
      gyro_bias_local_ids.pop();
    }

    local_window.addVariables(*graph, variable_ids);

    // we fix the first variables to remove degrees of freedom
    // especially important with slim factors
    local_window.variable(pose_local_ids.front())->setStatus(VariableBase::Status::Fixed);
    local_window.variable(vel_local_ids.front())->setStatus(VariableBase::Status::Fixed);
    local_window.variable(acc_bias_local_ids.front())->setStatus(VariableBase::Status::Fixed);
    local_window.variable(gyro_bias_local_ids.front())->setStatus(VariableBase::Status::Fixed);

    solver.setGraph(local_window);
    solver.compute();

    std::cout << solver.iterationStats() << std::endl;

    std::cout << "curr vel after opt: " << curr_vel_var->estimate().transpose() << "\n";
    std::cout << "curr bias_acc after opt: " << curr_bias_acc->estimate().transpose() << "\n";
    std::cout << "curr bias_gyro after opt: " << curr_bias_gyro->estimate().transpose() << "\n";

    prev_pose_var  = curr_pose_var;
    prev_vel_var   = curr_vel_var;
    prev_bias_acc  = curr_bias_acc;
    prev_bias_gyro = curr_bias_gyro;
  }
  // unfix all the variables and do final optimization
  for (auto it = graph->variables().begin(); it != graph->variables().end(); ++it) {
    VariableBase* v = it.value();
    if (v->graphId() > 4) {
      v->setStatus(VariableBase::Active);
    }
  }
  solver.setGraph(graph);
  solver.compute();

  std::cout << "final optimization" << std::endl;
  std::cout << solver.iterationStats() << std::endl;
  const std::string boss_graph_filename = example_folder + "/imu_gps_optimized.boss";

  solver.saveGraph(boss_graph_filename);
  std::cout << "graph written in " << boss_graph_filename << std::endl;

  ofstream dumper;
  dumper.open(example_folder + "/output.txt");

  for (size_t i = 0; i < graph->variables().size(); ++i) {
    const VarPoseImuType* curr_pose_var = dynamic_cast<const VarPoseImuType*>(graph->variable(i));
    if (curr_pose_var == nullptr) {
      continue;
    }

    const auto& curr_pose            = curr_pose_var->estimate();
    const Vector3f& pose_translation = curr_pose.translation();
    dumper << pose_translation.x() << " " << pose_translation.y() << " " << pose_translation.z()
           << " "
           << "\n";
  }
  dumper.close();

  cov_dump.close();
  det_dump.close();
}

void loadKittiData(KittiCalibration& kitti_calibration,
                   vector<ImuMeasurement>& imu_measurements,
                   vector<GpsMeasurement>& gps_measurements) {
  string line;

  // Read IMU metadata and compute relative sensor pose transforms
  // BodyPtx BodyPty BodyPtz BodyPrx BodyPry BodyPrz AccelerometerSigma GyroscopeSigma
  // IntegrationSigma AccelerometerBiasSigma GyroscopeBiasSigma AverageDelta
  ifstream imu_metadata(example_folder + "/data/KittiImuBiasedMetadata.txt");
  printf("-- Reading sensor metadata\n");

  getline(imu_metadata, line, '\n'); // ignore the first line
  // Load Kitti calibration
  getline(imu_metadata, line, '\n');
  sscanf(line.c_str(),
         "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
         &kitti_calibration.body_ptx,
         &kitti_calibration.body_pty,
         &kitti_calibration.body_ptz,
         &kitti_calibration.body_prx,
         &kitti_calibration.body_pry,
         &kitti_calibration.body_prz,
         &kitti_calibration.accelerometer_sigma,
         &kitti_calibration.gyroscope_sigma,
         &kitti_calibration.integration_sigma,
         &kitti_calibration.accelerometer_bias_sigma,
         &kitti_calibration.gyroscope_bias_sigma,
         &kitti_calibration.average_delta_t);
  printf("IMU metadata: %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
         kitti_calibration.body_ptx,
         kitti_calibration.body_pty,
         kitti_calibration.body_ptz,
         kitti_calibration.body_prx,
         kitti_calibration.body_pry,
         kitti_calibration.body_prz,
         kitti_calibration.accelerometer_sigma,
         kitti_calibration.gyroscope_sigma,
         kitti_calibration.integration_sigma,
         kitti_calibration.accelerometer_bias_sigma,
         kitti_calibration.gyroscope_bias_sigma,
         kitti_calibration.average_delta_t);

  // Read IMU data
  // Time dt accelX accelY accelZ omegaX omegaY omegaZ
  printf("-- Reading IMU measurements from file\n");
  {
    ifstream imu_data(example_folder + "/data/KittiImuBiased.txt");
    getline(imu_data, line, '\n'); // ignore the first line

    double time = 0, dt = 0, acc_x = 0, acc_y = 0, acc_z = 0, gyro_x = 0, gyro_y = 0, gyro_z = 0;
    while (!imu_data.eof()) {
      getline(imu_data, line, '\n');
      sscanf(line.c_str(),
             "%lf %lf %lf %lf %lf %lf %lf %lf",
             &time,
             &dt,
             &acc_x,
             &acc_y,
             &acc_z,
             &gyro_x,
             &gyro_y,
             &gyro_z);

      ImuMeasurement measurement;
      measurement.time          = time;
      measurement.dt            = dt;
      measurement.accelerometer = Vector3f(acc_x, acc_y, acc_z);
      measurement.gyroscope     = Vector3f(gyro_x, gyro_y, gyro_z);
      imu_measurements.push_back(measurement);
    }
  }

  // Read GPS data
  // Time,X,Y,Z
  printf("-- Reading GPS measurements from file\n");
  {
    ifstream gps_data(example_folder + "/data/KittiGps.txt");
    getline(gps_data, line, '\n'); // ignore the first line

    double time = 0, gps_x = 0, gps_y = 0, gps_z = 0;
    while (!gps_data.eof()) {
      getline(gps_data, line, '\n');
      sscanf(line.c_str(), "%lf,%lf,%lf,%lf", &time, &gps_x, &gps_y, &gps_z);

      GpsMeasurement measurement;
      measurement.time     = time;
      measurement.position = Vector3f(gps_x, gps_y, gps_z);
      gps_measurements.push_back(measurement);
    }
  }
}