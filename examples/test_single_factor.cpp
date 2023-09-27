#include "common/common.h"

#include "variables_and_factors/imu_preintegration_factor.h"

int main() {
  using namespace srrg2_core;
  using namespace srrg2_solver;
  using namespace test_imu;
  // using TrajectoryType = SE3EightTrajectory;

  /*   Solver solver;
    solver.param_termination_criteria.setValue(nullptr);
    solver.param_max_iterations.pushBack(30);
    IterationAlgorithmBasePtr alg(new IterationAlgorithmGN);
    solver.param_algorithm.setValue(alg);
    FactorGraphPtr graph(new FactorGraph);

    float T    = 10;
    float freq = 300;

    std::shared_ptr<TrajectoryType> traj = std::make_shared<TrajectoryType>(T);
    FakeImu imu(traj, freq, 102030);

    std::vector<std::pair<ImuMeasurement, srrg2_core::Isometry3f>> data;
    imu.generateData(data);

    // create two variables and a preintegration factor between them
    // then fix the first one and see if the second is optimized
    // to see if ADErrorFactor() is implemented correctly
    VariableSE3QuaternionRightAD *pose_start, *pose_end;
    VariableVector3AD *vel_start, *vel_end;
    VariableVector3AD *bias_acc_start, *bias_acc_end;
    VariableVector3AD *bias_gyro_start, *bias_gyro_end;

    pose_start = new VariableSE3QuaternionRightAD();
    pose_start->setStatus(VariableBase::Fixed);
    vel_start = new VariableVector3AD();
    vel_start->setStatus(VariableBase::Fixed);
    bias_acc_start = new VariableVector3AD();
    bias_acc_start->setStatus(VariableBase::Fixed);
    bias_gyro_start = new VariableVector3AD();
    // bias_gyro_start->setStatus(VariableBase::Fixed);

    pose_end      = new VariableSE3QuaternionRightAD();
    vel_end       = new VariableVector3AD();
    bias_acc_end  = new VariableVector3AD();
    bias_gyro_end = new VariableVector3AD();

    graph->addVariable(VariableBasePtr(pose_start));
    // graph->addVariable(VariableBasePtr(pose_end));
    graph->addVariable(VariableBasePtr(vel_start));
    // graph->addVariable(VariableBasePtr(vel_end));
    graph->addVariable(VariableBasePtr(bias_acc_start));
    graph->addVariable(VariableBasePtr(bias_gyro_start));
    // graph->addVariable(VariableBasePtr(bias_acc_end));
    // graph->addVariable(VariableBasePtr(bias_gyro_end)); */

  ImuPreintegrationFactorAD* imu_factor = new ImuPreintegrationFactorAD();

  // graph->addFactor(FactorBasePtr(imu_factor));

  /*  solver.setGraph(graph);
   solver.compute();

   // imu preintegration
   ImuPreintegrator integrator;
   integrator.reset();
   float dt = 1 / imu.freq();
   std::cout << "dt: " << dt << std::endl;

   float dT = 3.0;

   Isometry3f initial_pose;
   for (size_t i = 0; i < data.size(); ++i) {
     if (i == 0) {
       initial_pose        = data.at(i).second;
       ImuMeasurement meas = data.at(i).first;
       integrator.preintegrate(meas, dt);
       srrg2_core::Vector3f pos, vel, acc;
       imu.trajectory().sampleTrajectory(meas.timestamp, pos, vel, acc);
       pose_start->setEstimate(initial_pose);
       vel_start->setEstimate(vel);
       bias_acc_start->setEstimate(Vector3f::Zero());
       bias_gyro_start->setEstimate(Vector3f::Zero());
       continue;
     }
     ImuMeasurement meas = data.at(i).first;
     if (meas.timestamp > dT)
       break;
     integrator.preintegrate(meas, dt);
   } */

  return 0;
}