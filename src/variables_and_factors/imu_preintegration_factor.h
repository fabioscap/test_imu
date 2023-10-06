#pragma once

#include "common/common.h"

#include "imu_preintegrator/imu_preintegrator.h"
#include "imu_preintegrator/imu_preintegrator_slim.h"

#include <srrg_solver/solver_core/ad_error_factor.h>
#include <srrg_solver/solver_core/error_factor.h>

#include <srrg_solver/variables_and_factors/types_3d/variable_se3.h>
#include <srrg_solver/variables_and_factors/types_3d/variable_se3_ad.h>
#include <srrg_solver/variables_and_factors/types_common/all_types.h>

namespace srrg2_solver {

  class ImuPreintegrationFactorBase {
  public:
    // inline Measurement& measurement() {
    //   return pm_;
    //  }

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
    using dMatrix3f = srrg2_core::MatrixN_<srrg2_core::ad::DualValuef, 3>;
    using dVector3f = srrg2_core::Vector_<srrg2_core::ad::DualValuef, 3>;

    using ADErrorVectorType = typename BaseType::ADErrorVectorType;
    using VariableTupleType = typename BaseType::VariableTupleType;

    ADErrorVectorType operator()(VariableTupleType& vars);

    void setMeasurement(const test_imu::ImuPreintegrator& preintegrator);

    void _drawImpl(ViewerCanvasPtr canvas_) const override;

    inline void grav(const Vector3f& grav) {
      convertMatrix(grav_, grav);
    }

    // protected:
    dMatrix3f delta_R_;
    dVector3f delta_v_;
    dVector3f delta_p_;

    // nominal value for the biases
    dVector3f bias_acc_nom_;
    dVector3f bias_gyro_nom_;

    // bias correction matrices
    dMatrix3f dR_db_gyro_;
    dMatrix3f dv_db_acc_;
    dMatrix3f dv_db_gyro_;
    dMatrix3f dp_db_acc_;
    dMatrix3f dp_db_gyro_;

    DualValuef dT_;

    dVector3f grav_;
  };

  class ImuPreintegrationFactorSlimAD
    : public ADErrorFactor_<9,
                            VariableSE3QuaternionRightAD, // pose_from
                            VariableVector3AD,            // vel_from
                            VariableSE3QuaternionRightAD, // pose_to
                            VariableVector3AD>,           // vel_to

      public ImuPreintegrationFactorBase {
  public:
    // clang-format off
    using BaseType = ADErrorFactor_<9,
                                    VariableSE3QuaternionRightAD, // pose_from
                                    VariableVector3AD,            // vel_from
                                    VariableSE3QuaternionRightAD, // pose_to
                                    VariableVector3AD>;           // vel_to
    // clang-format on
    using dMatrix3f = srrg2_core::MatrixN_<srrg2_core::ad::DualValuef, 3>;
    using dVector3f = srrg2_core::Vector_<srrg2_core::ad::DualValuef, 3>;

    using ADErrorVectorType = typename BaseType::ADErrorVectorType;
    using VariableTupleType = typename BaseType::VariableTupleType;

    ADErrorVectorType operator()(VariableTupleType& vars);

    void setMeasurement(const test_imu::ImuPreintegratorSlim& preintegrator);

    inline void grav(const Vector3f& grav) {
      convertMatrix(grav_, grav);
    }

    // protected:
    dMatrix3f delta_R_;
    dVector3f delta_v_;
    dVector3f delta_p_;

    DualValuef dT_;

    dVector3f grav_;
  };
} // namespace srrg2_solver