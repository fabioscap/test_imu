#include <fstream>
#include <iostream>

#define _USE_MATH_DEFINES
#include "variables_and_factors/gps_factor.h"
#include "variables_and_factors/imu_preintegration_factor.h"
#include "variables_and_factors/instances.h"

#include "imu_preintegrator/imu_preintegrator.h"
#include "imu_preintegrator/imu_preintegrator_ukf.h"

#include "synthetic/synthetic.h"

#include "srrg_solver/solver_core/factor_graph.h"
#include "srrg_solver/solver_core/instances.h"

#include "srrg_solver/solver_core/internals/linear_solvers/instances.h"
#include "srrg_solver/solver_core/solver.h"
#include "srrg_solver/variables_and_factors/types_3d/all_types.h"
#include "srrg_solver/variables_and_factors/types_3d/instances.h"
#include <srrg_solver/solver_core/iteration_algorithm_dl.h>
#include <srrg_solver/solver_core/iteration_algorithm_lm.h>

#include "srrg_solver/solver_core/factor_graph_view.h"

#include <Eigen/SparseCore>

#include <thread>
#include <unistd.h>

#include <string.h>

#include <chrono>

#include <random>

using namespace std;
using namespace srrg2_core;
using namespace srrg2_solver;

const std::string example_folder("/workspace/src/test_imu");

using test_imu::ImuMeasurement;

struct Sigmas {
  float acc       = 0.00175f;
  float gyro      = 0.00175f;
  float bias_acc  = 0.000167f;
  float bias_gyro = 0.000167f;

  Vector3f gps = Vector3f(0.0002f, 0.0002f, 0.0002f);
} sigmas;

struct GpsMeasurement {
  double time;
  Vector3f position; // x,y,z
};

void generateData(std::vector<GpsMeasurement>& gps_measurements,
                  std::vector<test_imu::ImuMeasurement>& imu_measurements,
                  Vector3f&,
                  bool dump_gt);

bool string_in_array(const string& query, int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i)
    if (std::string(argv[i]) == query)
      return true;

  return false;
}

int main(int argc, char* argv[]) {
  // incremental optimization parameters
  const int window_size = 64;
  const int gps_skip    = 1;

  const int first_gps = 0;

  ofstream prediction_error("prediction_error.txt");

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

  if (string_in_array("ukf_slim", argc, argv)) {
    use_slim = true;
    use_ukf  = true;
  }

  // inspect covariances
  ofstream cov_dump("/workspace/src/test_imu/covariance.txt");
  // ofstream det_dump("/workspace/src/test_imu/determinant.txt");

  variables_and_factors_3d_registerTypes();
  variables_and_factors_imu_registerTypes();

  vector<ImuMeasurement> imu_measurements;
  vector<GpsMeasurement> gps_measurements;

  Vector3f vel0;
  generateData(gps_measurements, imu_measurements, vel0, true);

  using AlgorithmType = IterationAlgorithmLM;

  Solver solver;
  solver.param_termination_criteria.setValue(nullptr);
  solver.param_max_iterations.pushBack(50);
  IterationAlgorithmBasePtr alg(new AlgorithmType);
  std::shared_ptr<AlgorithmType> temp = std::dynamic_pointer_cast<AlgorithmType>(alg);

  // temp->param_damping.setValue(1000);
  //  temp->param_user_lambda_init.setValue(20);
  //  temp->param_variable_damping.setValue(false);

  solver.param_algorithm.setValue(alg);
  // solver.param_verbose.setValue(true);
  FactorGraphPtr graph(new FactorGraph);

  using VarPoseImuType    = VariableSE3ExpMapRight;
  using VarVelImuType     = VariableVector3;
  using ImuBiasVar        = VariableVector3;
  using FactorGpsType     = GpsFactor;
  using FactorBiasType    = BiasErrorFactor;
  using FactorImuType     = ImuPreintegrationFactor;
  using FactorImuSlimType = ImuPreintegrationFactorSlim;

  // initialization
  // hash map for variable timestamps
  std::vector<std::pair<double, size_t>> variable_timestamps;

  Isometry3f imu_in_body; // imu in gps_rtk
  imu_in_body.setIdentity();

  const Vector3f& init_gps_pose = gps_measurements.at(first_gps).position;
  Isometry3f init_pose          = Isometry3f::Identity();
  init_pose.translation()       = init_gps_pose;
  init_pose.linear()            = test_imu::Rz<float>(M_PI / 2);

  VarPoseImuType* prev_pose_var = new VarPoseImuType();
  prev_pose_var->setEstimate(init_pose);
  prev_pose_var->setGraphId(0);
  graph->addVariable(VariableBasePtr(prev_pose_var));
  variable_timestamps.push_back(
    std::make_pair(gps_measurements.at(first_gps).time, prev_pose_var->graphId()));

  FactorGpsType* gps_factor = new FactorGpsType();
  gps_factor->setVariableId(0, 0);

  gps_factor->setInformationMatrix(sigmas.gps.asDiagonal().inverse());
  gps_factor->setMeasurement(init_gps_pose);
  graph->addFactor(FactorBasePtr(gps_factor));

  // std::cout << "adding pose variable with id: " << prev_pose_var->graphId() << "\n";
  pose_local_ids.push(prev_pose_var->graphId());
  variable_ids.insert(prev_pose_var->graphId());

  VarVelImuType* prev_vel_var = new VarVelImuType();
  prev_vel_var->setEstimate(vel0);
  prev_vel_var->setGraphId(1);
  // prev_vel_var->setStatus(VariableBase::Status::Fixed);
  graph->addVariable(VariableBasePtr(prev_vel_var));
  vel_local_ids.push(prev_vel_var->graphId());
  variable_ids.insert(prev_vel_var->graphId());
  // std::cout << "adding vel variable with id: " << prev_vel_var->graphId() << "\n";

  ImuBiasVar* prev_bias_acc = new ImuBiasVar();
  prev_bias_acc->setGraphId(2);
  prev_bias_acc->setEstimate(0 * Vector3f::Ones());
  // prev_bias_acc->setStatus(VariableBase::Fixed);
  graph->addVariable(VariableBasePtr(prev_bias_acc));
  acc_bias_local_ids.push(prev_bias_acc->graphId());
  variable_ids.insert(prev_bias_acc->graphId());

  ImuBiasVar* prev_bias_gyro = new ImuBiasVar();
  prev_bias_gyro->setGraphId(3);
  prev_bias_gyro->setEstimate(0 * Vector3f::Ones());
  // prev_bias_gyro->setStatus(VariableBase::Fixed);
  graph->addVariable(VariableBasePtr(prev_bias_gyro));
  gyro_bias_local_ids.push(prev_bias_gyro->graphId());
  variable_ids.insert(prev_bias_gyro->graphId());

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

  imu_preintegrator->setNoiseGyroscope(Vector3f::Constant(sigmas.gyro));
  imu_preintegrator->setNoiseAccelerometer(Vector3f::Constant(sigmas.acc));
  imu_preintegrator->setNoiseBiasGyroscope(Vector3f::Constant(sigmas.bias_gyro));
  imu_preintegrator->setNoiseBiasAccelerometer(Vector3f::Constant(sigmas.bias_acc));

  size_t j                               = 0;
  size_t included_imu_measurements_count = 0;
  for (size_t i = first_gps + 1; i < gps_measurements.size(); i = i + gps_skip) {
    int prev_idx = (i == first_gps + 1) ? first_gps : i - gps_skip;

    const double t_previous = gps_measurements.at(prev_idx).time;
    const double gps_time   = gps_measurements.at(i).time;
    const double dT         = gps_time - t_previous;

    std::cout << std::fixed << "gps_idx: " << i << " t: " << gps_time << "\n";

    imu_preintegrator->reset();
    imu_preintegrator->setBiasAcc(prev_bias_acc->estimate());
    imu_preintegrator->setBiasGyro(prev_bias_gyro->estimate());

    // time to preintegrate
    auto start = std::chrono::high_resolution_clock::now();
    while (j < imu_measurements.size() && imu_measurements.at(j).timestamp <= gps_time) {
      if (imu_measurements.at(j).timestamp >= t_previous) {
        imu_preintegrator->preintegrate(imu_measurements.at(j),
                                        imu_measurements.at(j + 1).timestamp -
                                          imu_measurements.at(j).timestamp);
        included_imu_measurements_count++;
      }

      j++;
    }
    auto stop     = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

    std::cout << "Time taken: " << duration.count() << " nanoseconds" << std::endl;

    // order of var indeces
    // 1 pose from imu
    // 2 vel from imu
    // 3 pose to imu
    // 4 vel to imu

    if (imu_preintegrator->measurements().size() < 3)
      continue;

    const Vector3f curr_gps_pose = gps_measurements.at(i).position;
    Isometry3f curr_pose         = Isometry3f::Identity();
    curr_pose.translation()      = curr_gps_pose;
    curr_pose.linear()           = prev_pose_var->estimate().linear();

    float pred_err =
      (curr_gps_pose - prev_pose_var->estimate().translation() + prev_vel_var->estimate() * dT)
        .norm();
    // compare gps position with constant velocity prediction
    prediction_error << pred_err << "\n";

    VarPoseImuType* curr_pose_var = new VarPoseImuType();
    curr_pose_var->setEstimate(curr_pose);
    curr_pose_var->setGraphId(graph_id++);
    graph->addVariable(VariableBasePtr(curr_pose_var));
    variable_timestamps.push_back(
      std::make_pair(gps_measurements.at(i).time, curr_pose_var->graphId()));
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
    curr_bias_acc->setEstimate(prev_bias_acc->estimate());
    curr_bias_gyro->setEstimate(prev_bias_gyro->estimate());
    graph->addVariable(VariableBasePtr(curr_bias_acc));
    //  curr_bias_acc->setStatus(VariableBase::Fixed);
    acc_bias_local_ids.push(curr_bias_acc->graphId());
    variable_ids.insert(curr_bias_acc->graphId());

    graph->addVariable(VariableBasePtr(curr_bias_gyro));
    gyro_bias_local_ids.push(curr_bias_gyro->graphId());
    variable_ids.insert(curr_bias_gyro->graphId());
    //  curr_bias_gyro->setStatus(VariableBase::Fixed);
    std::cout << "imu preintegrator has absorbed: " << imu_preintegrator->measurements().size()
              << " measurements.\n";
    if (use_imu) {
      // std::cout << imu_preintegrator->sigma() << "\n\n";
      if (!use_slim) {
        FactorImuType* imu_factor = new FactorImuType();
        imu_factor->setOffset(imu_in_body);
        imu_factor->setGrav(Vector3f(0.f, 0.f, 0.0f));
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
        FactorImuSlimType* imu_factor = new FactorImuSlimType();
        imu_factor->setOffset(imu_in_body);
        imu_factor->setGrav(Vector3f(0.f, 0.f, 0.0f));
        imu_factor->setVariableId(0, prev_pose_var->graphId());
        imu_factor->setVariableId(1, prev_vel_var->graphId());
        imu_factor->setVariableId(2, curr_pose_var->graphId());
        imu_factor->setVariableId(3, curr_vel_var->graphId());
        imu_factor->setVariableId(4, prev_bias_acc->graphId());
        imu_factor->setVariableId(5, prev_bias_gyro->graphId());

        FactorBiasType* bias_factor = new FactorBiasType();
        bias_factor->setVariableId(0, prev_bias_acc->graphId());
        bias_factor->setVariableId(1, prev_bias_gyro->graphId());

        bias_factor->setVariableId(2, curr_bias_acc->graphId());
        bias_factor->setVariableId(3, curr_bias_gyro->graphId());
        Matrix6f bias_sigma = Matrix6f::Identity();
        bias_sigma.block<3, 3>(0, 0) *= dT * sigmas.bias_acc;
        bias_sigma.block<3, 3>(3, 3) *= dT * sigmas.bias_gyro;
        bias_factor->setInformationMatrix(bias_sigma.inverse());
        imu_factor->setMeasurement(*imu_preintegrator);

        graph->addFactor(FactorBasePtr(bias_factor));
        graph->addFactor(FactorBasePtr(imu_factor));
      }
    }
    FactorGpsType* gps_factor = new FactorGpsType();

    gps_factor->setVariableId(0, curr_pose_var->graphId());

    gps_factor->setInformationMatrix(sigmas.gps.asDiagonal().inverse());
    gps_factor->setMeasurement(curr_gps_pose);

    if (use_gps) {
      graph->addFactor(FactorBasePtr(gps_factor));
    }
    cov_dump << "cov:\n" << imu_preintegrator->sigma() << "\n\n";
    //  cov_dump << "omega:\n" << imu_preintegrator->sigma().inverse() << "\n\n";
    // det_dump << imu_preintegrator->sigma().determinant() << "\n";
    imu_preintegrator->reset();

    // local optimization
    FactorGraphView local_window;

    if (window_size > 0) {
      while (pose_local_ids.size() > window_size) {
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

      if (use_imu) {
        solver.setGraph(local_window);
        solver.compute();
      }
    } else {
      if (use_imu) {
        solver.setGraph(graph);
        solver.compute();
      }
    }

    // std::cout << solver.iterationStats() << std::endl;

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

  const std::string boss_graph_filename = "/workspace/src/test_imu/imu_gps_optimized.boss";

  if (!use_imu) {
    graph->setSerializationLevel(-1);
    graph->write(boss_graph_filename);
  } else {
    solver.setGraph(graph);
    solver.compute();
    std::cout << "final optimization" << std::endl;
    std::cout << solver.iterationStats() << std::endl;

    solver.saveGraph(boss_graph_filename);
    std::cout << "graph written in " << boss_graph_filename << std::endl;
  }
  // dumper.close();
  cov_dump.close();
  // det_dump.close();

  std::ofstream fp_out(std::string(PROJECT_FOLDER) + "/synthetic.csv");
  if (!fp_out)
    throw std::runtime_error("cannot open output csv file");

  for (size_t i = 0; i < variable_timestamps.size(); ++i) {
    auto pair      = variable_timestamps.at(i);
    double time    = pair.first;
    size_t pose_id = pair.second;

    VariableBase* var        = graph->variables().find(pose_id).value();
    VarPoseImuType* pose_var = dynamic_cast<VarPoseImuType*>(var);
    if (!pose_var) {
      throw std::runtime_error("error fetching variable");
    }

    Eigen::Isometry3f pose = pose_var->estimate();
    auto rpy               = pose.rotation().eulerAngles(0, 1, 2);

    fp_out << std::fixed << time << ", " << pose.translation().x() << ", " << pose.translation().y()
           << ", " << pose.translation().z() << ", " << rpy.x() << ", " << rpy.y() << ", "
           << rpy.z() << "\n";
  }

  fp_out.close();
}

void generateData(std::vector<GpsMeasurement>& gps_measurements,
                  std::vector<test_imu::ImuMeasurement>& imu_measurements,
                  Vector3f& vel_zero,
                  bool dump_gt) {
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

  std::vector<std::tuple<test_imu::ImuMeasurement, Vector3f, Isometry3f>> data;
  imu.generateData(data, true);

  std::random_device rd;
  std::mt19937 gen(rd());

  std::normal_distribution<double> distr(0, 1);

  imu_measurements.clear();
  gps_measurements.clear();

  std::ofstream pose_gt;
  if (dump_gt)
    pose_gt.open(std::string(PROJECT_FOLDER) + "/synthetic_gt.csv");

  for (size_t i = 0; i < data.size(); ++i) {
    // add the imu measurements to the vector
    imu_measurements.push_back(std::get<0>(data.at(i)));

    if (i % gps_freq == 0) {
      // add a new gps measurement
      Isometry3f& pose = std::get<2>(data.at(i));
      GpsMeasurement gps_meas;
      // get the time from IMU
      gps_meas.time = std::get<0>(data.at(i)).timestamp;
      // get translation and add noise
      gps_meas.position = pose.translation() + Vector3f(distr(gen) * sigmas.gps.x(),
                                                        distr(gen) * sigmas.gps.y(),
                                                        distr(gen) * sigmas.gps.z());
      gps_measurements.push_back(gps_meas);
      auto rpy = pose.rotation().eulerAngles(0, 1, 2);

      pose_gt << std::fixed << std::get<0>(data.at(i)).timestamp << ", " << pose.translation().x()
              << ", " << pose.translation().y() << ", " << pose.translation().z() << ", " << rpy.x()
              << ", " << rpy.y() << ", " << rpy.z() << "\n";
    }

    if (i == 0)
      vel_zero = std::get<1>(data.at(i));
  }
}