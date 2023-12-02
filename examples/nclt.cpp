#include <fstream>
#include <iostream>

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
#include <srrg_solver/solver_core/iteration_algorithm_dl.h>
#include <srrg_solver/solver_core/iteration_algorithm_gn.h>
#include <srrg_solver/solver_core/iteration_algorithm_lm.h>

#include "srrg_solver/solver_core/factor_graph_view.h"

using namespace std;
using namespace srrg2_core;
using namespace srrg2_solver;
using test_imu::ImuMeasurement;

const std::string data_folder("/workspace/src/test_imu/data/nclt");

struct GpsMeasurement {
  double time;
  Vector3f position; // x,y,z
};

struct Sigmas {
  float acc       = 0.001f;
  float gyro      = 0.000175f;
  float bias_acc  = 0.000167f;
  float bias_gyro = 2.91e-006f;

  Vector3f gps = Vector3f(1, 1, 9);
} sigmas;

bool string_in_array(const string& query, int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i)
    if (std::string(argv[i]) == query)
      return true;

  return false;
}

void loadNCLTData(vector<ImuMeasurement>& imu_measurements,
                  vector<GpsMeasurement>& gps_measurements);
int main(int argc, char* argv[]) {
  // incremental optimization parameters
  const int window_size = 16;
  const int gps_skip    = 4;

  const int first_gps = 4;

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

  // inspect covariances
  // ofstream cov_dump("/workspace/src/test_imu/covariance.txt");
  // ofstream det_dump("/workspace/src/test_imu/determinant.txt");

  variables_and_factors_3d_registerTypes();
  variables_and_factors_imu_registerTypes();

  vector<ImuMeasurement> imu_measurements;
  vector<GpsMeasurement> gps_measurements;
  loadNCLTData(imu_measurements, gps_measurements);

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
  // hash map for variable timestamps
  std::vector<std::pair<double, size_t>> variable_timestamps;

  Isometry3f imu_in_body; // imu in gps_rtk
  imu_in_body.setIdentity();
  imu_in_body.translation() << -0.11, -0.18, -0.71;

  const Vector3f& init_gps_pose = gps_measurements.at(first_gps).position;
  Isometry3f init_pose          = Isometry3f::Identity();
  init_pose.translation()       = init_gps_pose;
  // init_pose.linear()            = test_imu::Rx<float>(M_PI);

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
  prev_vel_var->setEstimate(Vector3f::Zero());
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

  size_t j = 0;
  for (size_t i = first_gps + 1; i < gps_measurements.size() / 8; i = i + gps_skip) {
    std::cout << "gps_idx: " << i << "\n";
    int prev_idx = (i == first_gps + 1) ? first_gps : i - gps_skip;

    const double t_previous = gps_measurements.at(prev_idx).time;
    const double gps_time   = gps_measurements.at(i).time;
    const double dT         = gps_time - t_previous;

    imu_preintegrator->reset();
    imu_preintegrator->setBiasAcc(prev_bias_acc->estimate());
    imu_preintegrator->setBiasGyro(prev_bias_gyro->estimate());

    std::cout << "integrating IMU measurements in [" << t_previous << "," << gps_time << "]\n";
    while (j < imu_measurements.size() && imu_measurements.at(j).timestamp <= gps_time) {
      if (imu_measurements.at(j).timestamp >= t_previous) {
        imu_preintegrator->preintegrate(imu_measurements.at(j),
                                        imu_measurements.at(j + 1).timestamp -
                                          imu_measurements.at(j).timestamp);
      }

      j++;
    }

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
    if (use_imu) {
      if (!use_slim) {
        ImuPreintegrationFactorAD* imu_factor = new ImuPreintegrationFactorAD();
        imu_factor->setOffset(imu_in_body);
        imu_factor->grav(Vector3f(0.f, 0.f, 9.80655));
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
        ImuPreintegrationFactorSlimAD* imu_factor = new ImuPreintegrationFactorSlimAD();
        imu_factor->setOffset(imu_in_body);
        imu_factor->grav(Vector3f(0.f, 0.f, 9.80655));
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
        bias_sigma.block<3, 3>(0, 0) *= dT * sigmas.bias_acc;
        bias_sigma.block<3, 3>(3, 3) *= dT * sigmas.bias_gyro;
        bias_factor->setInformationMatrix(bias_sigma.inverse());
        imu_factor->setMeasurement(*imu_preintegrator);

        graph->addFactor(FactorBasePtr(bias_factor));
        graph->addFactor(FactorBasePtr(imu_factor));
      }
    }
    FactorGpsType* gps_factor = new FactorGpsType();
    RobustifierBase* clamp    = new RobustifierClamp();
    clamp->param_chi_threshold.setValue(1000.0);
    // gps_factor->setRobustifier(clamp);

    gps_factor->setVariableId(0, curr_pose_var->graphId());

    gps_factor->setInformationMatrix(sigmas.gps.asDiagonal().inverse());
    gps_factor->setMeasurement(curr_gps_pose);

    if (use_gps) {
      graph->addFactor(FactorBasePtr(gps_factor));
    }
    // cov_dump << "cov:\n" << imu_preintegrator->sigma() << "\n\n";
    //  cov_dump << "omega:\n" << imu_preintegrator->sigma().inverse() << "\n\n";
    // det_dump << imu_preintegrator->sigma().determinant() << "\n";
    imu_preintegrator->reset();

    // local optimization
    FactorGraphView local_window;

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
  // cov_dump.close();
  // det_dump.close();

  std::ofstream fp_out(std::string(PROJECT_FOLDER) + "/nclt.csv");
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

    fp_out << time << ", " << pose.translation().x() << ", " << pose.translation().y() << ", "
           << pose.translation().z() << ", " << rpy.x() << ", " << rpy.y() << ", " << rpy.z()
           << "\n";
  }

  fp_out.close();
}

void loadNCLTData(vector<ImuMeasurement>& imu_measurements,
                  vector<GpsMeasurement>& gps_measurements) {
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
           "%lf,%*f,%*f,%*f,%f,%f,%f,%f,%f,%f",
           &meas.timestamp,
           &meas.acceleration.x(),
           &meas.acceleration.y(),
           &meas.acceleration.z(),
           &meas.angular_vel.x(),
           &meas.angular_vel.y(),
           &meas.angular_vel.z());
    meas.timestamp *= 1e-6;
    imu_measurements.push_back(meas);
  }

  ifstream gps_file(data_folder + "/" + gps_filename);
  if (!gps_file) {
    throw std::runtime_error("loadNCLTData: cannot open GPS data file at" + data_folder + "/" +
                             gps_filename + ".");
  }

  Eigen::Vector3f body_offset = Eigen::Vector3f::Zero();
  if (gps_filename == "gt.csv") {
    body_offset << +0.24, 0, +1.24;
  } else if (gps_filename == "gps.csv") {
    body_offset << 0, +0.25, +0.49;
  }
  std::cout << "gps body offset: " << body_offset.transpose() << "\n";

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
    if (gps_filename != "gt.csv") {
      sscanf(
        line.c_str(), "%lf,%d,%*d,%lf,%lf,%lf,%*f,%*f", &meas.time, &fix_mode, &lat, &lng, &alt);

      if (fix_mode >= 3) {
        // convert using formulas in pdf
        meas.position.x() = sin(lat - lat0) * rns;
        meas.position.y() = sin(lng - lng0) * rew * cos(lat0);
        meas.position.z() = alt0 - alt;
      } else {
        continue;
      }
    } else {
      sscanf(line.c_str(),
             "%lf,%f,%f,%f,%*f,%*f,%*f",
             &meas.time,
             &meas.position.x(),
             &meas.position.y(),
             &meas.position.z());
    }
    meas.position += body_offset;
    meas.time *= 1e-6;
    gps_measurements.push_back(meas);
  }
}
