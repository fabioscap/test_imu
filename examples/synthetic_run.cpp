#include <boost/program_options.hpp>
#include <chrono>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/stat.h>

#include <cstring>

#include "synthetic/synthetic.h"

#include <Eigen/Dense>

// srrg
#include "imu_preintegrator/imu_preintegrator.h"
#include "imu_preintegrator/imu_preintegrator_ukf.h"
#include "variables_and_factors/gps_factor.h"
#include "variables_and_factors/imu_preintegration_factor.h"
#include "variables_and_factors/instances.h"

#include "srrg_solver/solver_core/factor_graph.h"
#include "srrg_solver/solver_core/instances.h"

#include "srrg_solver/solver_core/solver.h"
#include "srrg_solver/variables_and_factors/types_3d/all_types.h"
#include "srrg_solver/variables_and_factors/types_3d/instances.h"

#include <srrg_solver/solver_core/iteration_algorithm_dl.h>
#include <srrg_solver/solver_core/iteration_algorithm_gn.h>
#include <srrg_solver/solver_core/iteration_algorithm_lm.h>

#include "srrg_solver/solver_core/factor_graph_view.h"

using test_imu::ImuMeasurement;

struct Sigmas {
  float acc;
  float gyro;
  float bias_acc;
  float bias_gyro;

  float gps;
} sigmas;

struct GpsMeasurement {
  double time;
  Eigen::Vector3f position; // x,y,z
};

// trajectory parameters
float T;
float imu_freq;
int gps_step;
int seed = 102030;
std::string trajectory_type;

// generated synthetic data
std::vector<ImuMeasurement> imu_measurements;
std::vector<GpsMeasurement> gps_measurements;
Eigen::Vector3f vel0;

std::string results_dir;

void parseArguments(int argc, char* argv[]);

void generateData(std::vector<GpsMeasurement>& gps_measurements,
                  std::vector<test_imu::ImuMeasurement>& imu_measurements,
                  Eigen::Vector3f& vel_zero,
                  bool dump_gt);

void run_srrg();

int main(int argc, char* argv[]) {
  // generate the folder in which results are stored
  auto now     = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::tm* ptm = std::localtime(&now);
  std::stringstream folder_ss;
  folder_ss << std::put_time(ptm, "%Y-%m-%d-%H-%M-%S");
  results_dir = PROJECT_FOLDER + std::string("/simulations/") + folder_ss.str();
  mkdir(results_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  // parse the arguments from the configuration file
  parseArguments(argc, argv);

  // generate synthetic data
  generateData(gps_measurements, imu_measurements, vel0, true);

  // run srrg
  run_srrg();

  return 0;
}

void generateData(std::vector<GpsMeasurement>& gps_measurements,
                  std::vector<test_imu::ImuMeasurement>& imu_measurements,
                  Eigen::Vector3f& vel_zero,
                  bool dump_gt) {
  using namespace test_imu;

  std::shared_ptr<Trajectory> traj;
  if (trajectory_type == "eight")
    traj = std::make_shared<SE3EightTrajectory>(T);
  else if (trajectory_type == "spike") {
    traj = std::make_shared<MeasurementTrajectory>(T);
  } else {
    traj = std::make_shared<SE3SineTrajectory>(T);
  }

  test_imu::FakeImu imu(traj, imu_freq, seed);

  imu.std_acc()       = sigmas.acc;
  imu.std_gyro()      = sigmas.gyro;
  imu.std_bias_acc()  = sigmas.bias_acc;
  imu.std_bias_gyro() = sigmas.bias_gyro;

  std::vector<std::tuple<test_imu::ImuMeasurement, Eigen::Vector3f, Eigen::Isometry3f>> data;
  imu.generateData(data, true);

  std::random_device rd;
  std::mt19937 gen(rd());

  std::normal_distribution<double> distr(0, 1);

  imu_measurements.clear();
  gps_measurements.clear();

  std::ofstream pose_gt;

  if (dump_gt)
    pose_gt.open(results_dir + "/synthetic_gt.csv");

  for (size_t i = 0; i < data.size(); ++i) {
    // add the imu measurements to the vector
    imu_measurements.push_back(std::get<0>(data.at(i)));

    if (i % gps_step == 0) {
      // add a new gps measurement
      Eigen::Isometry3f& pose = std::get<2>(data.at(i));
      GpsMeasurement gps_meas;
      // get the time from IMU
      gps_meas.time = std::get<0>(data.at(i)).timestamp;
      // get translation and add noise
      gps_meas.position =
        pose.translation() +
        Eigen::Vector3f(distr(gen) * sigmas.gps, distr(gen) * sigmas.gps, distr(gen) * sigmas.gps);
      gps_measurements.push_back(gps_meas);
    }
    if (dump_gt) {
      Eigen::Isometry3f& pose = std::get<2>(data.at(i));
      auto rpy                = pose.rotation().eulerAngles(0, 1, 2);
      pose_gt << std::fixed << std::get<0>(data.at(i)).timestamp << ", " << pose.translation().x()
              << ", " << pose.translation().y() << ", " << pose.translation().z() << ", " << rpy.x()
              << ", " << rpy.y() << ", " << rpy.z() << "\n";
    }
    if (i == 0)
      vel_zero = std::get<1>(data.at(i));
  }
}

void parseArguments(int argc, char* argv[]) {
  namespace po = boost::program_options;

  po::options_description desc("Allowed options");
  // clang-format off
  desc.add_options()("config-file,c", po::value<std::string>()->required(), "Configuration file")(
    // synthetic trajectory parameters
    "imu_frequency", po::value<float>()->default_value(50.0), "IMU frequency parameter in Hz")(
    "gps_step", po::value<float>()->default_value(50.0), "Insert a GPS correction every [gps_step] imu measures")(
    "T", po::value<float>()->default_value(100.0), "Trajectory period in seconds")(
    "trajectory", po::value<std::string>()->default_value("eight"), "Trajectory type")(

    // noise parameters
    "acc_std", po::value<float>()->default_value(0.00175f), "Accelerometer noise standard deviation")(
    "gyro_std", po::value<float>()->default_value(0.00175f), "Gyroscope noise standard deviation")(
    "acc_bias_std", po::value<float>()->default_value(0.000167f), "Accelerometer bias noise standard deviation")(
    "gyro_bias_std", po::value<float>()->default_value(0.000167f), "Gyroscope bias noise standard deviation")(
    "gps_std", po::value<float>()->default_value(0.000167f), "GPS position standard deviation")

    ;
  // clang-format on
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    exit(1);
  }

  if (vm.count("config-file")) {
    std::string config_filename = vm["config-file"].as<std::string>();
    std::ifstream config_file(config_filename);
    if (config_file.is_open()) {
      try {
        // parse the settings
        po::store(po::parse_config_file(config_file, desc), vm);
        po::notify(vm);
        // copy the settings used to run this simulation
        std::filesystem::copy(config_filename,
                              results_dir + std::string("/") +
                                std::filesystem::path(config_filename).filename().string());
      } catch (const std::exception& e) {
        std::cerr << "Error parsing config file: " << e.what() << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "Unable to open config file: " << config_filename << std::endl;
      exit(1);
    }
  }

  // set sigmas
  sigmas.acc       = vm["acc_std"].as<float>();
  sigmas.gyro      = vm["gyro_std"].as<float>();
  sigmas.bias_acc  = vm["acc_bias_std"].as<float>();
  sigmas.bias_gyro = vm["gyro_bias_std"].as<float>();
  sigmas.gps       = vm["gps_std"].as<float>();

  // trajectory parameters
  T = vm["T"].as<float>();
  // imu measurements /sec
  imu_freq = vm["imu_frequency"].as<float>();
  // insert a gps measurement every gps_step imu_measurements
  gps_step = vm["gps_step"].as<float>();
  // trajectory type
  trajectory_type = vm["trajectory"].as<std::string>();
}

void run_srrg() {
  using namespace test_imu;
  using namespace srrg2_core;
  using namespace srrg2_solver;

  variables_and_factors_3d_registerTypes();
  variables_and_factors_imu_registerTypes();

  std::ofstream log_file(results_dir + std::string("/log.txt"));

  using AlgorithmType = IterationAlgorithmLM;

  bool use_gps  = true;
  bool use_slim = true;
  bool use_ukf  = false;

  size_t window_size = 0;
  size_t first_gps   = 0;
  size_t gps_skip    = 1;

  Solver solver;
  solver.param_termination_criteria.setValue(nullptr);
  solver.param_max_iterations.pushBack(50);
  IterationAlgorithmBasePtr alg(new AlgorithmType);
  std::shared_ptr<AlgorithmType> temp = std::dynamic_pointer_cast<AlgorithmType>(alg);

  solver.param_algorithm.setValue(alg);
  FactorGraphPtr graph(new FactorGraph);

  // initialization
  // hash map for variable timestamps
  std::vector<std::pair<double, size_t>> variable_timestamps;

  IdSet variable_ids;
  std::queue<VariableBase::Id> pose_local_ids;
  std::queue<VariableBase::Id> vel_local_ids;
  std::queue<VariableBase::Id> acc_bias_local_ids;
  std::queue<VariableBase::Id> gyro_bias_local_ids;

  Isometry3f imu_in_body;
  imu_in_body.setIdentity();

  const Vector3f& init_gps_pose = gps_measurements.at(first_gps).position;
  Isometry3f init_pose          = Isometry3f::Identity();
  init_pose.translation()       = init_gps_pose;
  // init_pose.linear()            = test_imu::Rz<float>(M_PI / 2);

  VariableSE3ExpMapRight* prev_pose_var = new VariableSE3ExpMapRight();
  prev_pose_var->setEstimate(init_pose);
  prev_pose_var->setGraphId(0);
  graph->addVariable(VariableBasePtr(prev_pose_var));
  variable_timestamps.push_back(
    std::make_pair(gps_measurements.at(first_gps).time, prev_pose_var->graphId()));

  GpsFactor* gps_factor = new GpsFactor();
  gps_factor->setVariableId(0, 0);

  Matrix3f gps_information = Vector3f::Constant(sigmas.gps).asDiagonal().inverse();
  gps_factor->setInformationMatrix(gps_information);
  gps_factor->setMeasurement(init_gps_pose);
  graph->addFactor(FactorBasePtr(gps_factor));

  pose_local_ids.push(prev_pose_var->graphId());
  variable_ids.insert(prev_pose_var->graphId());

  VariableVector3* prev_vel_var = new VariableVector3();
  prev_vel_var->setEstimate(vel0);
  prev_vel_var->setGraphId(1);
  graph->addVariable(VariableBasePtr(prev_vel_var));
  vel_local_ids.push(prev_vel_var->graphId());
  variable_ids.insert(prev_vel_var->graphId());

  VariableVector3* prev_bias_acc = new VariableVector3();
  prev_bias_acc->setGraphId(2);
  prev_bias_acc->setEstimate(0 * Vector3f::Ones());
  graph->addVariable(VariableBasePtr(prev_bias_acc));
  acc_bias_local_ids.push(prev_bias_acc->graphId());
  variable_ids.insert(prev_bias_acc->graphId());

  VariableVector3* prev_bias_gyro = new VariableVector3();
  prev_bias_gyro->setGraphId(3);
  prev_bias_gyro->setEstimate(0 * Vector3f::Ones());
  graph->addVariable(VariableBasePtr(prev_bias_gyro));
  gyro_bias_local_ids.push(prev_bias_gyro->graphId());
  variable_ids.insert(prev_bias_gyro->graphId());

  size_t graph_id = 4;

  test_imu::ImuPreintegratorBase* imu_preintegrator;
  if (use_slim) {
    log_file << "slim ";
    if (use_ukf) {
      log_file << "ukf\n";
      imu_preintegrator = new test_imu::ImuPreintegratorUKFSlim();
    } else {
      log_file << "normal\n";
      imu_preintegrator = new test_imu::ImuPreintegratorSlim();
    }
  } else {
    log_file << "full ";
    if (use_ukf) {
      log_file << "ukf\n";
      imu_preintegrator = new test_imu::ImuPreintegratorUKF();
    } else {
      log_file << "normal\n";
      imu_preintegrator = new test_imu::ImuPreintegrator();
    }
  }

  imu_preintegrator->setNoiseGyroscope(Vector3f::Constant(sigmas.gyro * sigmas.gyro));
  imu_preintegrator->setNoiseAccelerometer(Vector3f::Constant(sigmas.acc * sigmas.acc));
  imu_preintegrator->setNoiseBiasGyroscope(Vector3f::Constant(sigmas.bias_gyro * sigmas.bias_gyro));
  imu_preintegrator->setNoiseBiasAccelerometer(
    Vector3f::Constant(sigmas.bias_acc * sigmas.bias_acc));

  size_t j                               = 0;
  size_t included_imu_measurements_count = 0;
  for (size_t i = first_gps + 1; i < gps_measurements.size(); i = i + gps_skip) {
    int prev_idx = (i == first_gps + 1) ? first_gps : i - gps_skip;

    const double t_previous = gps_measurements.at(prev_idx).time;
    const double gps_time   = gps_measurements.at(i).time;
    const double dT         = gps_time - t_previous;

    log_file << std::fixed << "gps_idx: " << i << " t: " << gps_time << "\n";

    imu_preintegrator->reset();
    imu_preintegrator->setBiasAcc(prev_bias_acc->estimate());
    imu_preintegrator->setBiasGyro(prev_bias_gyro->estimate());

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

    log_file << "Time taken: " << duration.count() << " nanoseconds" << std::endl;

    const Vector3f curr_gps_pose = gps_measurements.at(i).position;
    Isometry3f curr_pose         = Isometry3f::Identity();
    curr_pose.translation()      = curr_gps_pose;
    curr_pose.linear()           = prev_pose_var->estimate().linear();

    VariableSE3ExpMapRight* curr_pose_var = new VariableSE3ExpMapRight();
    curr_pose_var->setEstimate(curr_pose);
    curr_pose_var->setGraphId(graph_id++);
    graph->addVariable(VariableBasePtr(curr_pose_var));
    variable_timestamps.push_back(
      std::make_pair(gps_measurements.at(i).time, curr_pose_var->graphId()));
    pose_local_ids.push(curr_pose_var->graphId());
    variable_ids.insert(curr_pose_var->graphId());

    VariableVector3* curr_vel_var = new VariableVector3();
    curr_vel_var->setGraphId(graph_id++);
    const Vector3f prev_vel =
      static_cast<VariableVector3*>(graph->variable(graph_id - 3))->estimate();
    curr_vel_var->setEstimate(prev_vel);
    graph->addVariable(VariableBasePtr(curr_vel_var));
    vel_local_ids.push(curr_vel_var->graphId());
    variable_ids.insert(curr_vel_var->graphId());

    VariableVector3* curr_bias_acc  = new VariableVector3();
    VariableVector3* curr_bias_gyro = new VariableVector3();
    curr_bias_acc->setGraphId(graph_id++);
    curr_bias_gyro->setGraphId(graph_id++);
    curr_bias_acc->setEstimate(prev_bias_acc->estimate());
    curr_bias_gyro->setEstimate(prev_bias_gyro->estimate());

    graph->addVariable(VariableBasePtr(curr_bias_acc));
    acc_bias_local_ids.push(curr_bias_acc->graphId());
    variable_ids.insert(curr_bias_acc->graphId());

    graph->addVariable(VariableBasePtr(curr_bias_gyro));
    gyro_bias_local_ids.push(curr_bias_gyro->graphId());
    variable_ids.insert(curr_bias_gyro->graphId());

    if (!use_slim) {
      ImuPreintegrationFactor* imu_factor = new ImuPreintegrationFactor();
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
      ImuPreintegrationFactorSlim* imu_factor = new ImuPreintegrationFactorSlim();
      imu_factor->setOffset(imu_in_body);
      imu_factor->setGrav(Vector3f(0.f, 0.f, 0.0f));
      imu_factor->setVariableId(0, prev_pose_var->graphId());
      imu_factor->setVariableId(1, prev_vel_var->graphId());
      imu_factor->setVariableId(2, curr_pose_var->graphId());
      imu_factor->setVariableId(3, curr_vel_var->graphId());
      imu_factor->setVariableId(4, prev_bias_acc->graphId());
      imu_factor->setVariableId(5, prev_bias_gyro->graphId());

      BiasErrorFactor* bias_factor = new BiasErrorFactor();
      bias_factor->setVariableId(0, prev_bias_acc->graphId());
      bias_factor->setVariableId(1, prev_bias_gyro->graphId());

      bias_factor->setVariableId(2, curr_bias_acc->graphId());
      bias_factor->setVariableId(3, curr_bias_gyro->graphId());
      Matrix6f bias_sigma = Matrix6f::Identity();
      bias_sigma.block<3, 3>(0, 0) *= sigmas.bias_acc / dT;
      bias_sigma.block<3, 3>(3, 3) *= dT * sigmas.bias_gyro / dT;
      bias_factor->setInformationMatrix(bias_sigma.inverse());
      imu_factor->setMeasurement(*imu_preintegrator);

      graph->addFactor(FactorBasePtr(bias_factor));
      graph->addFactor(FactorBasePtr(imu_factor));
    }

    GpsFactor* gps_factor = new GpsFactor();

    gps_factor->setVariableId(0, curr_pose_var->graphId());

    gps_factor->setInformationMatrix(gps_information);
    gps_factor->setMeasurement(curr_gps_pose);

    if (use_gps) {
      graph->addFactor(FactorBasePtr(gps_factor));
    }

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

      solver.setGraph(local_window);
      solver.compute();

    } else {
      solver.setGraph(graph);
      solver.compute();
    }

    log_file << solver.iterationStats() << std::endl;

    log_file << "curr vel after opt: " << curr_vel_var->estimate().transpose() << "\n";
    log_file << "curr bias_acc after opt: " << curr_bias_acc->estimate().transpose() << "\n";
    log_file << "curr bias_gyro after opt: " << curr_bias_gyro->estimate().transpose() << "\n";

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

  const std::string boss_graph_filename = results_dir + std::string("/graph.boss");

  solver.setGraph(graph);
  solver.compute();
  log_file << "final optimization" << std::endl;
  log_file << solver.iterationStats() << std::endl;

  solver.saveGraph(boss_graph_filename);
  log_file << "graph written in " << boss_graph_filename << std::endl;

  std::ofstream fp_out(results_dir + "/synthetic.csv");
  if (!fp_out)
    throw std::runtime_error("cannot open output csv file");

  for (size_t i = 0; i < variable_timestamps.size(); ++i) {
    auto pair      = variable_timestamps.at(i);
    double time    = pair.first;
    size_t pose_id = pair.second;

    VariableBase* var                = graph->variables().find(pose_id).value();
    VariableSE3ExpMapRight* pose_var = dynamic_cast<VariableSE3ExpMapRight*>(var);
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
