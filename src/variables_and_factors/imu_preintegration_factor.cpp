#include "imu_preintegration_factor.h"

#include <srrg_solver/solver_core/ad_error_factor_impl.cpp>
#include <srrg_solver/solver_core/error_factor_impl.cpp>

#include <srrg_solver/solver_core/instance_macros.h>

#define P_i 0
#define PHI_i 3
#define V_i 6
#define P_j 9
#define PHI_j 12
#define V_j 15
#define BA_i 18
#define BG_i 21
#define BA_j 24
#define BG_j 27

namespace srrg2_solver {

  template class ImuPreintegrationFactorBase<float,
                                             ErrorFactor_<15,
                                                          VariableSE3ExpMapRight, // pose_from
                                                          VariableVector3,        // vel_from
                                                          VariableSE3ExpMapRight, // pose_to
                                                          VariableVector3,        // vel_to
                                                          VariableVector3,        // bias_acc_from
                                                          VariableVector3,        // bias_gyro_from
                                                          VariableVector3,        // bias_acc_to
                                                          VariableVector3>>;
  template class ImuPreintegrationFactorBase<ad::DualValuef,
                                             ADErrorFactor_<15,
                                                            VariableSE3ExpMapRightAD, // pose_from
                                                            VariableVector3AD,        // vel_from
                                                            VariableSE3ExpMapRightAD, // pose_to
                                                            VariableVector3AD,        // vel_to
                                                            VariableVector3AD, // bias_acc_from
                                                            VariableVector3AD, // bias_gyro_from
                                                            VariableVector3AD, // bias_acc_to
                                                            VariableVector3AD>>;
  template class ImuPreintegrationFactorBase<float,
                                             ErrorFactor_<9,
                                                          VariableSE3ExpMapRight, // pose_from
                                                          VariableVector3,        // vel_from
                                                          VariableSE3ExpMapRight, // pose_to
                                                          VariableVector3,        // vel_to
                                                          VariableVector3,        // bias_from
                                                          VariableVector3>>;
  template class ImuPreintegrationFactorBase<ad::DualValuef,
                                             ADErrorFactor_<9,
                                                            VariableSE3ExpMapRightAD, // pose_from
                                                            VariableVector3AD,        // vel_from
                                                            VariableSE3ExpMapRightAD, // pose_to
                                                            VariableVector3AD,        // vel_to
                                                            VariableVector3AD,        // bias_from
                                                            VariableVector3AD>>;

  template <typename Scalar_>
  Matrix3_<Scalar_> jacobianExpMapSO3inv(const Vector3_<Scalar_> omega_) {
    Matrix3_<Scalar_> Jri;
    Scalar_ theta_square = omega_.dot(omega_);
    Scalar_ theta        = sqrt(theta_square);
    Matrix3_<Scalar_> W  = geometry3d::skew(omega_);
    if (theta_square < Scalar_(1e-8)) {
      Jri = Matrix3_<Scalar_>::Identity() + Scalar_(0.5) * W;
    } else {
      Jri = Matrix3_<Scalar_>::Identity() + Scalar_(0.5) * W +
            ((1 / theta_square) + (1 + cos(theta)) / (2 * theta * sin(theta))) * W * W;
    }
    return Jri;
  }

  template <typename Scalar_, typename BaseType_>
  template <size_t i, bool ad>
  inline auto srrg2_solver::ImuPreintegrationFactorBase<Scalar_, BaseType_>::getEstimateAt(
    typename BaseType::VariableTupleType& vars) {
    if constexpr (ad) {
      return vars.template at<i>()->adEstimate();
    } else {
      return vars.template at<i>()->estimate();
    }
  }

  template <typename Scalar_, typename BaseType_>
  auto ImuPreintegrationFactorBase<Scalar_, BaseType_>::computeResiduals(
    typename BaseType_::VariableTupleType& vars,
    bool error_only_) {
    // find out which type of factor is calling this function
    constexpr bool slim = (BaseType_::ErrorDim == 9);
    constexpr bool ad   = std::is_same<Scalar_, DualValuef>::value;

    // get the variables
    const Isometry3 T_i = getEstimateAt<0, ad>(vars) * offset_;
    const Vector3 v_i =
      getEstimateAt<1, ad>(vars) + offset_.linear() * (omega_i_.cross(offset_.translation()));
    const Isometry3 T_j = getEstimateAt<2, ad>(vars) * offset_;
    const Vector3 v_j =
      getEstimateAt<3, ad>(vars) + offset_.linear() * (omega_j_.cross(offset_.translation()));
    const Vector3 ba_i = getEstimateAt<4, ad>(vars);
    const Vector3 bg_i = getEstimateAt<5, ad>(vars);

    const Matrix3& R_i = T_i.linear();
    const Vector3& p_i = T_i.translation();
    const Matrix3& R_j = T_j.linear();
    const Vector3& p_j = T_j.translation();

    Vector_<Scalar_, BaseType_::ErrorDim> error;

    // bias correction via first order approximation around nominal values
    const Vector3 delta_bacc  = ba_i - bias_acc_nom_;
    const Vector3 delta_bgyro = bg_i - bias_gyro_nom_;

    // orientation error
    error.template segment<3>(0) = geometry3d::logMapSO3(
      ((delta_R_ * geometry3d::expMapSO3((dR_db_gyro_ * delta_bgyro).eval())).transpose() *
       R_i.transpose() * R_j)
        .eval());
    // velocity error
    error.template segment<3>(3) = R_i.transpose() * (v_j - v_i - grav_ * dT_) -
                                   (delta_v_ + dv_db_acc_ * delta_bacc + dv_db_gyro_ * delta_bgyro);
    // position error
    error.template segment<3>(6) =
      R_i.transpose() * (p_j - p_i - v_i * dT_ - DualValuef(0.5) * grav_ * dT_ * dT_) -
      (delta_p_ + dp_db_acc_ * delta_bacc + dp_db_gyro_ * delta_bgyro);

    if constexpr (!slim) {
      // complete the error vector
      const Vector3 ba_j            = getEstimateAt<6, ad>(vars);
      const Vector3 bg_j            = getEstimateAt<7, ad>(vars);
      error.template segment<3>(9)  = ba_j - ba_i;
      error.template segment<3>(12) = bg_j - bg_i;
    }

    // AD factors stop here
    if constexpr (ad) {
      return error;
    } else {
      // side effect on _e
      BaseType::_e = error;

      // now we compute jacobians (if necessary)
      if (error_only_) {
        return error;
      }

      // check if this is necessary
      // std::cout <<jacobian<0>() << "\n";
      BaseType::_J.setZero();

      // jacobians of delta_p_
      auto Jpi  = BaseType::template jacobian<0>();
      auto Jvi  = BaseType::template jacobian<1>();
      auto Jpj  = BaseType::template jacobian<2>();
      auto Jvj  = BaseType::template jacobian<3>();
      auto Jbai = BaseType::template jacobian<4>();
      auto Jbgi = BaseType::template jacobian<5>();

      if constexpr (!slim) {
        BaseType::template jacobian<6>().setZero();
        BaseType::template jacobian<7>().setZero();
        auto Jbaj                        = BaseType::template jacobian<6>();
        auto Jbgj                        = BaseType::template jacobian<7>();
        Jbaj.template block<3, 3>(9, 0)  = Matrix3::Identity();
        Jbgj.template block<3, 3>(12, 0) = Matrix3::Identity();
        // now the jacobians related to the bias propagation
        Jbai.template block<3, 3>(9, 0)  = -Matrix3::Identity();
        Jbgi.template block<3, 3>(12, 0) = -Matrix3::Identity();
      }

      // jacobians related to delta_p_, which starts at the 6th row
      Jpi.template block<3, 3>(6, 3) = geometry3d::skew(
        (R_i.transpose() * (p_j - p_i - v_i * dT_ - 0.5 * grav_ * dT_ * dT_)).eval());

      Jpi.template block<3, 3>(6, 0) = -Matrix3::Identity();

      Jvi.template block<3, 3>(6, 0) = -R_i.transpose() * dT_;

      Jbai.template block<3, 3>(6, 0) = -dp_db_acc_;

      Jpj.template block<3, 3>(6, 0) = R_i.transpose() * R_j;

      Jbgi.template block<3, 3>(6, 0) = -dp_db_gyro_;

      // jacobians related to delta_v_
      Jpi.template block<3, 3>(3, 3) =
        geometry3d::skew((R_i.transpose() * (v_j - v_i - grav_ * dT_)).eval());

      Jvi.template block<3, 3>(3, 0) = -R_i.transpose();

      Jbai.template block<3, 3>(3, 0) = -dv_db_acc_;

      Jvj.template block<3, 3>(3, 0) = R_i.transpose();

      Jbgi.template block<3, 3>(3, 0) = -dv_db_gyro_;

      // jacobians related to delta_R_
      Vector3 r_phi                  = error.head(3);
      Matrix3 Jinv                   = jacobianExpMapSO3inv(r_phi);
      Jpi.template block<3, 3>(0, 3) = -Jinv * R_j.transpose() * R_i;

      Jpj.template block<3, 3>(0, 3) = Jinv;

      Jbgi.template block<3, 3>(0, 0) =
        -Jinv * geometry3d::expMapSO3(r_phi).transpose() *
        geometry3d::jacobianExpMapSO3((dR_db_gyro_ * delta_bgyro).eval()) * dR_db_gyro_;

      // offset jacobians
      const Matrix3& Rt  = offset_.linear().transpose();
      const Vector3& tbi = offset_.translation();

      J_pert_.setIdentity();
      J_pert_.template block<3, 3>(P_i, P_i)     = Rt;
      J_pert_.template block<3, 3>(P_i, PHI_i)   = -Rt * geometry3d::skew(tbi);
      J_pert_.template block<3, 3>(PHI_i, PHI_i) = Rt;
      J_pert_.template block<3, 3>(V_i, PHI_i) =
        -R_i * geometry3d::skew(omega_i_.cross(tbi).eval());

      J_pert_.template block<3, 3>(P_j, P_j)     = Rt;
      J_pert_.template block<3, 3>(P_j, PHI_j)   = -Rt * geometry3d::skew(tbi);
      J_pert_.template block<3, 3>(PHI_j, PHI_j) = Rt;
      J_pert_.template block<3, 3>(V_j, PHI_j) =
        -R_j * geometry3d::skew(omega_j_.cross(tbi).eval());

      BaseType::_J = BaseType::_J * J_pert_;
    }

    return error; // suppress warning
  }

  template <typename Scalar_, typename BaseType_>
  void srrg2_solver::ImuPreintegrationFactorBase<Scalar_, BaseType_>::setMeasurement(
    test_imu::ImuPreintegratorBase& preintegrator) {
    using srrg2_core::ad::convertMatrix;

    const test_imu::BiasJacobians& bias_J_ptr = preintegrator.biasJacobians();

    // TODO is it bad to call convertMatrix even if the type is the same?
    convertMatrix(delta_R_, preintegrator.delta_R());
    convertMatrix(delta_v_, preintegrator.delta_v());
    convertMatrix(delta_p_, preintegrator.delta_p());
    convertMatrix(bias_acc_nom_, preintegrator.bias_acc());
    convertMatrix(bias_gyro_nom_, preintegrator.bias_gyro());
    convertMatrix(dR_db_gyro_, bias_J_ptr.dR_db_gyro);
    convertMatrix(dv_db_acc_, bias_J_ptr.dv_db_acc);
    convertMatrix(dv_db_gyro_, bias_J_ptr.dv_db_gyro);
    convertMatrix(dp_db_acc_, bias_J_ptr.dp_db_acc);
    convertMatrix(dp_db_gyro_, bias_J_ptr.dp_db_gyro);

    dT_ = Scalar_(preintegrator.dT());
    BaseType_::setInformationMatrix(preintegrator.sigma().inverse());

    convertMatrix(
      omega_i_,
      (preintegrator.measurements().front().angular_vel - preintegrator.bias_gyro()).eval());
    convertMatrix(
      omega_j_,
      (preintegrator.measurements().back().angular_vel - preintegrator.bias_gyro()).eval());
  }

  template <typename Scalar_, typename BaseType_>
  void ImuPreintegrationFactorBase<Scalar_, BaseType_>::setOffset(const Isometry3f& offset) {
    srrg2_core::ad::convertMatrix(offset_, offset);

    /*     constexpr bool slim = (BaseType_::ErrorDim == 9);
        constexpr bool ad   = std::is_same<Scalar_, DualValuef>::value;
     */
  }

  void ImuPreintegrationFactor::errorAndJacobian(bool error_only_) {
    computeResiduals(_variables, error_only_);
  }

  void ImuPreintegrationFactorSlim::errorAndJacobian(bool error_only_) {
    computeResiduals(_variables, error_only_);
  }

  ImuPreintegrationFactorAD::ADErrorVectorType
  ImuPreintegrationFactorAD::operator()(ImuPreintegrationFactorAD::VariableTupleType& vars) {
    return computeResiduals(vars, true);
  }

  ImuPreintegrationFactorSlimAD::ADErrorVectorType ImuPreintegrationFactorSlimAD::operator()(
    ImuPreintegrationFactorSlimAD::VariableTupleType& vars) {
    return computeResiduals(vars, true);
  }

  void BiasErrorFactor::errorAndJacobian(bool error_only_) {
    const Vector3f& bias_acc_i  = _variables.at<0>()->estimate();
    const Vector3f& bias_gyro_i = _variables.at<1>()->estimate();
    const Vector3f& bias_acc_j  = _variables.at<2>()->estimate();
    const Vector3f& bias_gyro_j = _variables.at<3>()->estimate();

    _e.segment<3>(0) = bias_acc_j - bias_acc_i;
    _e.segment<3>(3) = bias_gyro_j - bias_gyro_i;

    if (error_only_)
      return;

    auto Jbai = jacobian<0>(); // 15x3
    auto Jbgi = jacobian<1>(); // 15x3
    auto Jbaj = jacobian<2>(); // 15x3
    auto Jbgj = jacobian<3>(); // 15x3

    // now the jacobians related to the bias propagation
    Jbaj.block<3, 3>(0, 0).setIdentity();
    Jbgj.block<3, 3>(3, 0).setIdentity();

    Jbai.block<3, 3>(0, 0) = -Matrix3f::Identity();
    Jbgi.block<3, 3>(3, 0) = -Matrix3f::Identity();
  }

  BiasErrorFactorAD::ADErrorVectorType BiasErrorFactorAD::operator()(VariableTupleType& vars) {
    const dVector3f& bias_acc_i  = vars.at<0>()->adEstimate();
    const dVector3f& bias_gyro_i = vars.at<1>()->adEstimate();
    const dVector3f& bias_acc_j  = vars.at<2>()->adEstimate();
    const dVector3f& bias_gyro_j = vars.at<3>()->adEstimate();

    ADErrorVectorType error;

    error.segment<3>(0) = bias_acc_j - bias_acc_i;
    error.segment<3>(3) = bias_gyro_j - bias_gyro_i;

    return error;
  }

  template <typename Scalar_, typename BaseType_>
  void ImuPreintegrationFactorBase<Scalar_, BaseType_>::_drawImpl(ViewerCanvasPtr canvas_) const {
    Vector3f coords[2];
    coords[0] = reinterpret_cast<const VariableSE3ExpMapRight*>(BaseType_::variable(0))
                  ->estimate()
                  .translation();
    coords[1] = reinterpret_cast<const VariableSE3ExpMapRight*>(BaseType_::variable(2))
                  ->estimate()
                  .translation();

    float lw = 0.5;
    if (fabs(BaseType_::variableId(0) - BaseType_::variableId(1)) == 1) {
      lw *= 2;
    }
    lw *= (BaseType_::level() * 3 + 1);
    canvas_->pushColor();
    canvas_->pushLineWidth();
    canvas_->setLineWidth(lw);
    float fading   = 1. - 0.5 * BaseType_::level();
    Vector3f color = srrg2_core::ColorPalette::color3fBlue() * fading;
    canvas_->setColor(color);
    canvas_->putLine(2, coords);
    canvas_->popAttribute();
    canvas_->popAttribute();
  }

  INSTANTIATE(ImuPreintegrationFactor)
  INSTANTIATE(ImuPreintegrationFactorAD)
  INSTANTIATE(ImuPreintegrationFactorSlim)
  INSTANTIATE(ImuPreintegrationFactorSlimAD)
  INSTANTIATE(BiasErrorFactor)
  INSTANTIATE(BiasErrorFactorAD)

} // namespace srrg2_solver