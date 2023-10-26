#include "common/common.h"

#include "srrg_solver/solver_core/factor.h"
#include "srrg_solver/solver_core/factor_graph.h"
#include "srrg_solver/solver_core/iteration_algorithm_dl.h"
#include "srrg_solver/solver_core/iteration_algorithm_gn.h"
#include "srrg_solver/solver_core/iteration_algorithm_lm.h"
#include "srrg_solver/solver_core/solver.h"

#include "synthetic/synthetic.h"
#include "variables_and_factors/imu_preintegration_factor.h"

#include "imu_preintegrator/imu_preintegrator.h"
#include "imu_preintegrator/imu_preintegrator_ukf.h"

#include <type_traits>
#include <typeinfo>

#include <fstream>

template <typename Scalar>
void dumpCov(const Eigen::Matrix<Scalar, -1, -1>& cov);

// create two variables and a preintegration factor between them
// then fix the first one and see if the second is optimized
// to see if ADErrorFactor() is implemented correctly
int main() {
  using namespace srrg2_core;
  using namespace srrg2_solver;
  using namespace test_imu;
  using TrajectoryType = SE3EightTrajectory;

  constexpr bool slim = true;
  using PreintegratorType =
    std::conditional<slim, ImuPreintegratorUKFSlim, ImuPreintegratorUKF>::type;
  using FactorType =
    std::conditional<slim, ImuPreintegrationFactorSlimAD, ImuPreintegrationFactorAD>::type;
  using VariablePoseType = VariableSE3ExpMapRightAD;

  Solver solver;
  solver.param_termination_criteria.setValue(nullptr);
  solver.param_max_iterations.pushBack(100);
  IterationAlgorithmDL* alg(new IterationAlgorithmDL);
  // alg->param_damping.setValue(1);
  solver.param_algorithm.setValue(IterationAlgorithmBasePtr(alg));
  FactorGraphPtr graph(new FactorGraph);

  float T    = 100;
  float freq = 500;

  std::shared_ptr<TrajectoryType> traj = std::make_shared<TrajectoryType>(T);
  FakeImu imu(traj, freq, 102030);

  std::vector<std::tuple<ImuMeasurement, core::Vector3f, core::Isometry3f>> data;
  imu.generateData(data, false);

  VariablePoseType *pose_start, *pose_end;
  VariableVector3AD *vel_start, *vel_end;
  VariableVector3AD *bias_acc_start, *bias_acc_end;
  VariableVector3AD *bias_gyro_start, *bias_gyro_end;

  pose_start = new VariablePoseType();
  vel_start  = new VariableVector3AD();
  if (!slim) {
    bias_acc_start  = new VariableVector3AD();
    bias_gyro_start = new VariableVector3AD();
  }

  pose_end = new VariablePoseType();
  vel_end  = new VariableVector3AD();

  if (!slim) {
    bias_acc_end  = new VariableVector3AD();
    bias_gyro_end = new VariableVector3AD();
  }

  pose_start->setStatus(VariableBase::Fixed);
  vel_start->setStatus(VariableBase::Fixed);
  if (!slim) {
    bias_acc_start->setStatus(VariableBase::Fixed);
    bias_gyro_start->setStatus(VariableBase::Fixed);
    // bias_acc_end->setStatus(VariableBase::Fixed);
    // bias_gyro_end->setStatus(VariableBase::Fixed);
  }
  graph->addVariable(VariableBasePtr(pose_start));
  graph->addVariable(VariableBasePtr(pose_end));

  graph->addVariable(VariableBasePtr(vel_start));
  graph->addVariable(VariableBasePtr(vel_end));
  if (!slim) {
    graph->addVariable(VariableBasePtr(bias_acc_start));
    graph->addVariable(VariableBasePtr(bias_gyro_start));
    graph->addVariable(VariableBasePtr(bias_acc_end));
    graph->addVariable(VariableBasePtr(bias_gyro_end));
  }

  FactorType* imu_factor = new FactorType();
  imu_factor->setGraphId(0);

  imu_factor->setVariableId(0, pose_start->graphId());
  imu_factor->setVariableId(1, vel_start->graphId());
  imu_factor->setVariableId(2, pose_end->graphId());
  imu_factor->setVariableId(3, vel_end->graphId());
  if (!slim) {
    imu_factor->setVariableId(4, bias_acc_start->graphId());
    imu_factor->setVariableId(5, bias_gyro_start->graphId());
    imu_factor->setVariableId(6, bias_acc_end->graphId());
    imu_factor->setVariableId(7, bias_gyro_end->graphId());
  }

  graph->addFactor(FactorBasePtr(imu_factor));

  solver.setGraph(graph);

  // imu preintegration
  PreintegratorType* integrator_ptr = new PreintegratorType();
  PreintegratorType& integrator     = *integrator_ptr;
  integrator.reset();
  float dt = 1 / imu.freq();
  std::cout << "dt: " << dt << std::endl;

  float dT = 10 * dt + dt / 2;

  Isometry3f initial_pose;
  ImuMeasurement meas;
  size_t i;
  for (i = 0; i < data.size() - 1; ++i) {
    if (i == 0) {
      initial_pose = std::get<2>(data.at(i));
      meas         = std::get<0>(data.at(i));
      integrator.preintegrate(meas, dt);
      pose_start->setEstimate(initial_pose);
      vel_start->setEstimate(std::get<1>(data.at(i)));
      if (!slim) {
        bias_acc_start->setEstimate(Vector3f::Zero());
        bias_gyro_start->setEstimate(Vector3f::Zero());
      }

      continue;
    }
    meas = std::get<0>(data.at(i));

    if (meas.timestamp > dT)
      break;
    integrator.preintegrate(meas, dt);
  }

  imu_factor->setMeasurement(integrator);
  imu_factor->grav(Vector3f(0, 0, 0));

  solver.compute();
  std::cout << "after compute\n";
  std::cerr << solver.iterationStats() << std::endl;

  std::cout << "sigma:\n" << integrator.sigma() << "\n";
  std::cout << "sigma determinant: " << integrator.sigma().determinant() << "\n";
  std::cout << "omega: \n" << imu_factor->informationMatrix() << "\n";
  std::cout << "omega determinant: \n" << imu_factor->informationMatrix().determinant() << "\n";

  if (!slim) {
    std::cout << "sigma gay determinant:\n"
              << integrator.sigma().block<9, 9>(0, 0).determinant() << "\n";
  }
  std::cout << "gt:\n" << std::get<2>(data.at(i)).matrix() << "\n";

  std::cout << "preintegrated measurements:\n " << integrator.delta_R() << "\n"
            << integrator.delta_p() << "\n";

  std::cout << "estimate:\n" << pose_end->estimate().matrix() << "\n";

  std::cout << "rel pose t2v:\n"
            << geometry3d::t2v(Isometry3f(std::get<2>(data.at(i)).matrix().inverse() *
                                          pose_end->estimate().matrix()))
            << "\n";

  dumpCov(integrator.sigma());

  return 0;
}

template <typename Scalar>
void dumpCov(const Eigen::Matrix<Scalar, -1, -1>& cov) {
  // Save the covariance matrix to a data file
  std::ofstream datafile("/workspace/src/test_imu/examples/covariance_matrix.txt");
  datafile << cov;
  datafile.close();
}