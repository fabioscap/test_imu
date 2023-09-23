#include "common/common.h"

#include "imu_preintegrator/imu_preintegrator.h"

#include <srrg_solver/solver_core/ad_error_factor.h>
#include <srrg_solver/solver_core/error_factor.h>

#include <srrg_solver/variables_and_factors/types_3d/variable_se3.h>
#include <srrg_solver/variables_and_factors/types_3d/variable_se3_ad.h>
#include <srrg_solver/variables_and_factors/types_common/all_types.h>

namespace srrg2_solver {

  class ImuPreintegrationFactorBase {
    using Measurement = test_imu::PreintegratedImuMeasurement;

  public:
    // inline Measurement& measurement() {
    //   return pm_;
    //  }

  protected:
    // Measurement pm_;

    // gravity vector
    core::Vector3f grav_;
  };

  class ImuPreintegrationFactor : public ErrorFactor_<15,
                                                      VariableSE3QuaternionRight, // pose_from
                                                      VariableVector3,            // vel_from
                                                      VariableSE3QuaternionRight, // pose_to
                                                      VariableVector3,            // vel_to
                                                      VariableVector3,            // bias_acc_from
                                                      VariableVector3,            // bias_gyro_from
                                                      VariableVector3,            // bias_acc_to
                                                      VariableVector3>,           // bias_gyro_to,
                                  public ImuPreintegrationFactorBase {
  public:
    // get variable: _variables.at<idx>()
    // store error in _e
    // store jacobian de/dvidx in jacobian<idx>()

    // see if I can subclass pose-pose constraint

    void errorAndJacobian(bool error_only) override;

  protected:
  };

  class ImuPreintegrationFactorAD : public ADErrorFactor_<15,
                                                          VariableSE3QuaternionRightAD, // pose_from
                                                          VariableVector3AD,            // vel_from
                                                          VariableSE3QuaternionRightAD, // pose_to
                                                          VariableVector3AD,            // vel_to
                                                          VariableVector3AD,  // bias_acc_from
                                                          VariableVector3AD,  // bias_acc_to
                                                          VariableVector3AD,  // bias_gyro_from
                                                          VariableVector3AD>, // bias_gyro_to,
                                    public ImuPreintegrationFactorBase {
  public:
    // clang-format off
    using BaseType = ADErrorFactor_<15,
                                    VariableSE3QuaternionRightAD, // pose_from
                                    VariableVector3AD,            // vel_from
                                    VariableSE3QuaternionRightAD, // pose_to
                                    VariableVector3AD,            // vel_to
                                    VariableVector3AD,            // bias_acc_from
                                    VariableVector3AD,            // bias_acc_to
                                    VariableVector3AD,            // bias_gyro_from
                                    VariableVector3AD>;           // bias_gyro_to
    // clang-format on
    using ADErrorVectorType = BaseType::ADErrorVectorType;
    using VariableTupleType = BaseType::VariableTupleType;

    ADErrorVectorType operator()(VariableTupleType& vars) override;

  protected:
    // TODO
    // preintegrated measurement must become dual we do not want to copy from integrator to
    // measurement AND call convertMatrix
    // 1) blast PrientegratedImuMeasurement and create
    // setMeasurement(const ImuPreintegrator&)
    //
    // 2) embed PreintegratedImuMeasurement into integrator and make it do side effect on that
    // (maybe better)
  };

} // namespace srrg2_solver