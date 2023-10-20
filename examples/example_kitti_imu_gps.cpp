#include <iostream>
#define _USE_MATH_DEFINES
#include "variables_and_factors/imu_preintegration_factor.h"
#include "variables_and_factors/instances.h"

#include "srrg_solver/solver_core/factor_graph.h"
#include "srrg_solver/solver_core/instances.h"

#include "srrg_solver/solver_core/internals/linear_solvers/instances.h"
#include "srrg_solver/solver_core/internals/sparse_block_matrix/matrix_block_factory.h"
#include "srrg_solver/solver_core/solver.h"
#include "srrg_solver/variables_and_factors/types_3d/all_types.h"
#include "srrg_solver/variables_and_factors/types_3d/instances.h"
#include <srrg_solver/solver_core/iteration_algorithm_dl.h>
#include <srrg_system_utils/parse_command_line.h>
#include <srrg_system_utils/shell_colors.h>

#include <srrg_qgl_viewport/viewer_core_shared_qgl.h>
#include <srrg_solver/solver_core/factor_graph.h>

#include <thread>
#include <unistd.h>

#include <string.h>

using namespace std;
using namespace srrg2_core;
using namespace srrg2_solver;

using namespace srrg2_qgl_viewport;

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

int main(int argc, char* argv[]) {
  bool use_ukf = false;
  bool use_imu = true;
  if (argc > 1) {
    if (!std::strcmp(argv[1], "ukf"))
      use_ukf = true;
    else if (!std::strcmp(argv[1], "noimu"))
      use_imu = false;
  }
  variables_and_factors_3d_registerTypes();
  variables_and_factors_imu_registerTypes();

  KittiCalibration kitti_calibration;
  vector<ImuMeasurement> imu_measurements;
  vector<GpsMeasurement> gps_measurements;
  loadKittiData(kitti_calibration, imu_measurements, gps_measurements);

  Solver solver;
  // solver.param_termination_criteria.setValue(nullptr);
  solver.param_max_iterations.pushBack(15);
  IterationAlgorithmBasePtr alg(new IterationAlgorithmDL);
  solver.param_algorithm.setValue(alg);
  FactorGraphPtr graph(new FactorGraph);

  using VarPoseImuType = VariableSE3QuaternionRightAD;
  using VarVelImuType  = VariableVector3AD;
  using ImuBiasVar     = VariableVector3AD;
  using FactorImuType  = ImuPreintegrationFactorAD;
  using FactorGpsType  = GpsErrorFactorAD;

  // initialization

  const Vector3f& init_gps_pose = gps_measurements.at(0).position;
  Isometry3f init_pose          = Isometry3f::Identity();
  init_pose.translation()       = init_gps_pose;

  VarPoseImuType* prev_pose_var = new VarPoseImuType();
  prev_pose_var->setEstimate(init_pose);
  prev_pose_var->setGraphId(0);

  // instead of fixing the first variable, we put a gps factor (translation prior)
  // prev_pose_var->setStatus(VariableBase::Status::Fixed);
  FactorGpsType* gps_factor = new FactorGpsType();
  gps_factor->setVariableId(0, 0);
  float sigma_gps = 0.2;
  float info_gps  = 1 / (sigma_gps * sigma_gps);
  gps_factor->setInformationMatrix(Matrix3f::Identity() * info_gps);
  gps_factor->setMeasurement(init_gps_pose);

  graph->addVariable(VariableBasePtr(prev_pose_var));

  VarVelImuType* prev_vel_var = new VarVelImuType();
  prev_vel_var->setEstimate(Vector3f::Zero());
  prev_vel_var->setGraphId(1);
  prev_vel_var->setStatus(VariableBase::Status::Fixed);
  graph->addVariable(VariableBasePtr(prev_vel_var));

  ImuBiasVar* prev_bias_acc = new ImuBiasVar();
  prev_bias_acc->setGraphId(2);
  prev_bias_acc->setEstimate(Vector3f::Zero());
  prev_bias_acc->setStatus(VariableBase::Fixed);
  graph->addVariable(VariableBasePtr(prev_bias_acc));

  ImuBiasVar* prev_bias_gyro = new ImuBiasVar();
  prev_bias_gyro->setGraphId(3);
  prev_bias_gyro->setEstimate(Vector3f::Zero());
  prev_bias_gyro->setStatus(VariableBase::Fixed);
  graph->addVariable(VariableBasePtr(prev_bias_gyro));

  size_t graph_id = 4;

  test_imu::ImuPreintegratorBase* imu_preintegrator;
  if (use_ukf)
    imu_preintegrator = new test_imu::ImuPreintegratorUKF();
  else
    imu_preintegrator = new test_imu::ImuPreintegrator();

  imu_preintegrator->setNoiseGyroscope(Vector3f::Constant(kitti_calibration.gyroscope_sigma));
  imu_preintegrator->setNoiseAccelerometer(
    Vector3f::Constant(kitti_calibration.accelerometer_sigma));
  imu_preintegrator->setNoiseBiasGyroscope(
    Vector3f::Constant(kitti_calibration.gyroscope_bias_sigma));
  imu_preintegrator->setNoiseBiasAccelerometer(
    Vector3f::Constant(kitti_calibration.accelerometer_bias_sigma));

  size_t j = 0;
  for (size_t i = 1; i < gps_measurements.size(); ++i) {
    const double t_previous = gps_measurements.at(i - 1).time;
    const double gps_time   = gps_measurements.at(i).time;

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

    VarPoseImuType* curr_pose_var = new VarPoseImuType();
    curr_pose_var->setEstimate(curr_pose);
    curr_pose_var->setGraphId(graph_id++);
    graph->addVariable(VariableBasePtr(curr_pose_var));

    VarVelImuType* curr_vel_var = new VarVelImuType();
    curr_vel_var->setGraphId(graph_id++);
    const Vector3f prev_vel =
      static_cast<VarVelImuType*>(graph->variable(graph_id - 3))->estimate();
    curr_vel_var->setEstimate(prev_vel);
    graph->addVariable(VariableBasePtr(curr_vel_var));

    ImuBiasVar* curr_bias_acc  = new ImuBiasVar();
    ImuBiasVar* curr_bias_gyro = new ImuBiasVar();
    curr_bias_acc->setGraphId(graph_id++);
    curr_bias_gyro->setGraphId(graph_id++);
    curr_bias_acc->setEstimate(Vector3f::Zero());
    curr_bias_gyro->setEstimate(Vector3f::Zero());
    graph->addVariable(VariableBasePtr(curr_bias_acc));
    graph->addVariable(VariableBasePtr(curr_bias_gyro));

    FactorImuType* imu_factor = new FactorImuType();
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

    FactorGpsType* gps_factor = new FactorGpsType();
    gps_factor->setVariableId(0, curr_pose_var->graphId());
    float sigma_gps = 0.8;
    float info_gps  = 1 / (sigma_gps * sigma_gps);
    gps_factor->setInformationMatrix(Matrix3f::Identity() * info_gps);
    gps_factor->setMeasurement(curr_gps_pose);

    if (use_imu)
      graph->addFactor(FactorBasePtr(imu_factor));
    graph->addFactor(FactorBasePtr(gps_factor));
    imu_preintegrator->reset();

    if (i > 5) {
      solver.setGraph(graph);
      solver.compute();
    }
    std::cerr << solver.iterationStats() << std::endl;

    prev_pose_var  = curr_pose_var;
    prev_vel_var   = curr_vel_var;
    prev_bias_acc  = curr_bias_acc;
    prev_bias_gyro = curr_bias_gyro;

    // solver compute
  }

  solver.setGraph(graph);
  solver.compute();

  std::cerr << "final optimization" << std::endl;
  std::cerr << solver.iterationStats() << std::endl;
  const std::string boss_graph_filename = example_folder + "/imu_gps_optimized.boss";

  solver.saveGraph(boss_graph_filename);
  std::cerr << "graph written in " << boss_graph_filename << std::endl;

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

  /*   QApplication qapp(argc, argv);
    ViewerCoreSharedQGL viewer_core(argc, argv, &qapp);

    std::thread graph_t(viewGraph, viewer_core.getCanvas("viewer_core_shared_canvas"));
    viewer_core.startViewerServer();

    graph_t.join(); */
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

std::map<std::string, FactorGraphPtr> graphs;

bool gl_list_generated = false;
void viewGraph(ViewerCanvasPtr canvas_) {
  std::cerr << "canvas_ok" << std::endl;
  while (ViewerCoreSharedQGL::isRunning()) {
    if (!canvas_->_setup()) {
      usleep(10000);
      continue;
    }
    if (!gl_list_generated) {
      for (auto& g_it : graphs) {
        std::string name = g_it.first;
        auto graph       = g_it.second;
        canvas_->createList(name);
        canvas_->beginList(name);
        for (auto v : graph->variables()) {
          v.second->_drawImpl(canvas_);
        }

        for (auto f : graph->factors()) {
          f.second->_drawImpl(canvas_);
        }
        canvas_->endList();
      }
      gl_list_generated = true;
    } else {
      for (auto& g_it : graphs) {
        std::string name = g_it.first;
        canvas_->callList(name);
      }
    }
    canvas_->flush();
    sleep(1);
  }
}