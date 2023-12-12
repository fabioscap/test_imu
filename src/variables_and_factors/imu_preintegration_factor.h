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

  template <typename Scalar_, typename BaseType_>
  class ImuPreintegrationFactorBase : public BaseType_ {
  public:
    using BaseType = BaseType_;

    using Matrix3 = srrg2_core::MatrixN_<Scalar_, 3>;
    using Vector3 = srrg2_core::Vector_<Scalar_, 3>;
    using Vector9 = srrg2_core::Vector_<Scalar_, 9>;

    using Isometry3 = srrg2_core::Isometry3_<Scalar_>;

    void setMeasurement(test_imu::ImuPreintegratorBase& preintegrator);

    inline void setOffset(const Isometry3f& offset) {
      srrg2_core::ad::convertMatrix(offset_, offset);
    }

    inline void setGrav(const Vector3f& grav) {
      srrg2_core::ad::convertMatrix(grav_, grav);
    }

  protected:
    void _drawImpl(ViewerCanvasPtr canvas_) const override;
    void computeResiduals(const Isometry3& Ti,
                          const Vector3& vi,
                          const Isometry3& Tj,
                          const Vector3& vj,
                          const Vector3& ba,
                          const Vector3& bg) const;

    Matrix3 delta_R_;
    Vector3 delta_v_;
    Vector3 delta_p_;

    // nominal value for the biases
    Vector3 bias_acc_nom_;
    Vector3 bias_gyro_nom_;

    // bias correction matrices
    Matrix3 dR_db_gyro_;
    Matrix3 dv_db_acc_;
    Matrix3 dv_db_gyro_;
    Matrix3 dp_db_acc_;
    Matrix3 dp_db_gyro_;

    Scalar_ dT_;

    Vector3 grav_;

    Isometry3 offset_ = Isometry3::Identity(); // imu in body

    // correction jacobians
    // how a perturbation in body frame transforms into imu frame
    MatrixN_<Scalar_, 6> pert_J_ = MatrixN_<Scalar_, 6>::Zero();
  };

  class ImuPreintegrationFactor
    : public ImuPreintegrationFactorBase<float,
                                         ErrorFactor_<15,
                                                      VariableSE3ExpMapRight, // pose_from
                                                      VariableVector3,        // vel_from
                                                      VariableSE3ExpMapRight, // pose_to
                                                      VariableVector3,        // vel_to
                                                      VariableVector3,        // bias_acc_from
                                                      VariableVector3,        // bias_gyro_from
                                                      VariableVector3,        // bias_acc_to
                                                      VariableVector3>>       // bias_gyro_to,
  {
  public:
    void errorAndJacobian(bool error_only_ = false) override;
  };

  class ImuPreintegrationFactorAD
    : public ImuPreintegrationFactorBase<ad::DualValuef,
                                         ADErrorFactor_<15,
                                                        VariableSE3ExpMapRightAD, // pose_from
                                                        VariableVector3AD,        // vel_from
                                                        VariableSE3ExpMapRightAD, // pose_to
                                                        VariableVector3AD,        // vel_to
                                                        VariableVector3AD,        // bias_acc_from
                                                        VariableVector3AD,        // bias_gyro_from
                                                        VariableVector3AD,        // bias_acc_to
                                                        VariableVector3AD>>       // bias_gyro_to,
  {
  public:
    using ADErrorVectorType = typename BaseType::ADErrorVectorType;
    using VariableTupleType = typename BaseType::VariableTupleType;

    ADErrorVectorType operator()(VariableTupleType& vars);
  };

  class ImuPreintegrationFactorSlim
    : public ImuPreintegrationFactorBase<float,
                                         ErrorFactor_<9,
                                                      VariableSE3ExpMapRight, // pose_from
                                                      VariableVector3,        // vel_from
                                                      VariableSE3ExpMapRight, // pose_to
                                                      VariableVector3,        // vel_to
                                                      VariableVector3,        // bias_from
                                                      VariableVector3>>

  {
  public:
    void errorAndJacobian(bool error_only_ = false) override;
  };

  class ImuPreintegrationFactorSlimAD
    : public ImuPreintegrationFactorBase<ad::DualValuef,
                                         ADErrorFactor_<9,
                                                        VariableSE3ExpMapRightAD, // pose_from
                                                        VariableVector3AD,        // vel_from
                                                        VariableSE3ExpMapRightAD, // pose_to
                                                        VariableVector3AD,        // vel_to
                                                        VariableVector3AD,        // bias_from
                                                        VariableVector3AD>> {
  public:
    using ADErrorVectorType = typename BaseType::ADErrorVectorType;
    using VariableTupleType = typename BaseType::VariableTupleType;

    ADErrorVectorType operator()(VariableTupleType& vars);
  };

  class BiasErrorFactor : public ErrorFactor_<6,
                                              VariableVector3, // bias acc from
                                              VariableVector3, // bias_gyro_from
                                              VariableVector3, // bias acc to
                                              VariableVector3> // bias_gyro_to

  {
  public:
    using BaseType = ErrorFactor_<6,
                                  VariableVector3,  // bias acc from
                                  VariableVector3,  // bias_gyro_from
                                  VariableVector3,  // bias acc to
                                  VariableVector3>; // bias_gyro_to
    void errorAndJacobian(bool error_only_ = false) override;

    inline void _drawImpl(ViewerCanvasPtr canvas_) const override{};
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