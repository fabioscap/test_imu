#pragma once

#include "common/common.h"

#include "imu_preintegrator/imu_preintegrator_base.h"

#include <srrg_solver/solver_core/ad_error_factor.h>
#include <srrg_solver/solver_core/error_factor.h>

// #include <srrg_solver/variables_and_factors/types_3d/variable_se3.h>
// #include <srrg_solver/variables_and_factors/types_3d/variable_se3_ad.h>
#include <srrg_solver/variables_and_factors/types_common/all_types.h>

#include "variable_se3_expmap_right.h"

namespace srrg2_solver {

  class ImuPreintegrationFactorBase {
  public:
    // inline Measurement& measurement() {
    //   return pm_;
    //  }

  protected:
  };

  class ImuPreintegrationFactorAD : public ADErrorFactor_<15,
                                                          VariableSE3ExpMapRightAD, // pose_from
                                                          VariableVector3AD,        // vel_from
                                                          VariableSE3ExpMapRightAD, // pose_to
                                                          VariableVector3AD,        // vel_to
                                                          VariableVector3AD,        // bias_acc_from
                                                          VariableVector3AD,  // bias_gyro_from
                                                          VariableVector3AD,  // bias_acc_to
                                                          VariableVector3AD>, // bias_gyro_to,
                                    public ImuPreintegrationFactorBase {
  public:
    using BaseType = ADErrorFactor_<15,
                                    VariableSE3ExpMapRightAD, // pose_from
                                    VariableVector3AD,        // vel_from
                                    VariableSE3ExpMapRightAD, // pose_to
                                    VariableVector3AD,        // vel_to
                                    VariableVector3AD,        // bias_acc_from
                                    VariableVector3AD,        // bias_gyro_from
                                    VariableVector3AD,        // bias_acc_to
                                    VariableVector3AD>;       // bias_gyro_to

    using dMatrix3f   = srrg2_core::MatrixN_<srrg2_core::ad::DualValuef, 3>;
    using dVector3f   = srrg2_core::Vector_<srrg2_core::ad::DualValuef, 3>;
    using dIsometry3f = srrg2_core::Isometry3_<srrg2_core::ad::DualValuef>;

    using ADErrorVectorType = typename BaseType::ADErrorVectorType;
    using VariableTupleType = typename BaseType::VariableTupleType;

    ADErrorVectorType operator()(VariableTupleType& vars);

    void setMeasurement(test_imu::ImuPreintegratorBase& preintegrator);

    void _drawImpl(ViewerCanvasPtr canvas_) const override;

    inline void grav(const Vector3f& grav) {
      convertMatrix(grav_, grav);
    }

    inline void setOffset(const Isometry3f& offset) {
      // convertMatrix(offset_, offset);
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

    dIsometry3f offset_ = dIsometry3f::Identity(); // imu in body
  };

  class ImuPreintegrationFactorSlimAD : public ADErrorFactor_<9,
                                                              VariableSE3ExpMapRightAD, // pose_from
                                                              VariableVector3AD,        // vel_from
                                                              VariableSE3ExpMapRightAD, // pose_to
                                                              VariableVector3AD>,       // vel_to

                                        public ImuPreintegrationFactorBase {
  public:
    using BaseType = ADErrorFactor_<9,
                                    VariableSE3ExpMapRightAD, // pose_from
                                    VariableVector3AD,        // vel_from
                                    VariableSE3ExpMapRightAD, // pose_to
                                    VariableVector3AD>;       // vel_to

    using dMatrix3f   = srrg2_core::MatrixN_<srrg2_core::ad::DualValuef, 3>;
    using dVector3f   = srrg2_core::Vector_<srrg2_core::ad::DualValuef, 3>;
    using dIsometry3f = srrg2_core::Isometry3_<srrg2_core::ad::DualValuef>;

    using ADErrorVectorType = typename BaseType::ADErrorVectorType;
    using VariableTupleType = typename BaseType::VariableTupleType;

    ADErrorVectorType operator()(VariableTupleType& vars);

    void setMeasurement(test_imu::ImuPreintegratorBase& preintegrator);

    inline void grav(const Vector3f& grav) {
      convertMatrix(grav_, grav);
    }

    void _drawImpl(ViewerCanvasPtr canvas_) const override;

    inline void setOffset(const Isometry3f& offset) {
      convertMatrix(offset_, offset);
    }

    // protected:
    dMatrix3f delta_R_;
    dVector3f delta_v_;
    dVector3f delta_p_;

    DualValuef dT_;

    dVector3f grav_;

    dIsometry3f offset_ = dIsometry3f::Identity(); // imu in body
  };

  class BiasErrorFactorAD : public ADErrorFactor_<6,
                                                  VariableVector3AD, // bias acc from
                                                  VariableVector3AD, // bias_gyro_from
                                                  VariableVector3AD, // bias acc to
                                                  VariableVector3AD> // bias_gyro_to

  {
  public:
    using BaseType = ADErrorFactor_<6,
                                    VariableVector3AD,  // bias acc from
                                    VariableVector3AD,  // bias_gyro_from
                                    VariableVector3AD,  // bias acc to
                                    VariableVector3AD>; // bias_gyro_to

    using dMatrix3f = srrg2_core::MatrixN_<srrg2_core::ad::DualValuef, 3>;
    using dVector3f = srrg2_core::Vector_<srrg2_core::ad::DualValuef, 3>;

    using ADErrorVectorType = typename BaseType::ADErrorVectorType;
    using VariableTupleType = typename BaseType::VariableTupleType;

    ADErrorVectorType operator()(VariableTupleType& vars);

    inline void _drawImpl(ViewerCanvasPtr canvas_) const override{};
  };
} // namespace srrg2_solver