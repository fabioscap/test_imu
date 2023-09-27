#include "imu_preintegration_factor.h"

// include these two files otherwise nothing works and pia foco
#include <srrg_solver/solver_core/ad_error_factor_impl.cpp>
#include <srrg_solver/solver_core/error_factor_impl.cpp>

#include <srrg_solver/solver_core/instance_macros.h>

namespace srrg2_solver {

  ImuPreintegrationFactorAD::ADErrorVectorType
  ImuPreintegrationFactorAD::operator()(ImuPreintegrationFactorAD::VariableTupleType& vars) {
    const Isometry3_<DualValuef>& T_i = vars.at<0>()->adEstimate();
    const dMatrix3f& R_i              = T_i.linear();
    const dVector3f& p_i              = T_i.translation();
    const dVector3f& v_i              = vars.at<1>()->adEstimate();

    const Isometry3_<DualValuef>& T_j = vars.at<2>()->adEstimate();
    const dMatrix3f& R_j              = T_j.linear();
    const dVector3f& p_j              = T_j.translation();
    const dVector3f& v_j              = vars.at<3>()->adEstimate();

    const dVector3f& bias_acc_i  = vars.at<4>()->adEstimate();
    const dVector3f& bias_acc_j  = vars.at<5>()->adEstimate();
    const dVector3f& bias_gyro_i = vars.at<6>()->adEstimate();
    const dVector3f& bias_gyro_j = vars.at<7>()->adEstimate();

    ADErrorVectorType error;

    const dVector3f delta_bacc  = bias_acc_i - bias_acc_nom_;
    const dVector3f delta_bgyro = bias_gyro_i - bias_gyro_nom_;

    // we can also use any pose-pose error in place of the one in the paper

    // orientation error
    error.segment<3>(0) = geometry3d::logMapSO3(
      (delta_R_ * geometry3d::expMapSO3((dR_db_gyro_ * delta_bgyro).eval()) * R_i.transpose() * R_j)
        .eval());
    // velocity error
    error.segment<3>(3) = R_i.transpose() * (v_j - v_i - grav_ * dT_) -
                          (delta_v_ + dv_db_acc_ * delta_bacc + dv_db_gyro_ * delta_bgyro);
    // position error
    error.segment<3>(6) =
      R_i.transpose() * (p_j - p_i - v_i * dT_ - DualValuef(0.5) * grav_ * dT_ * dT_) -
      (delta_p_ + dp_db_acc_ * delta_bacc + dp_db_gyro_ * delta_bgyro);

    return error;
  }

  void ImuPreintegrationFactorAD::setMeasurement(const test_imu::ImuPreintegrator& preintegrator) {
    convertMatrix(delta_R_, preintegrator.delta_R());
    convertMatrix(delta_v_, preintegrator.delta_v());
    convertMatrix(delta_p_, preintegrator.delta_p());

    convertMatrix(bias_acc_nom_, preintegrator.bias_acc());
    convertMatrix(bias_gyro_nom_, preintegrator.bias_gyro());

    convertMatrix(dR_db_gyro_, preintegrator.dR_db_gyro());
    convertMatrix(dv_db_acc_, preintegrator.dv_db_acc());
    convertMatrix(dv_db_gyro_, preintegrator.dv_db_gyro());
    convertMatrix(dp_db_acc_, preintegrator.dp_db_acc());
    convertMatrix(dp_db_gyro_, preintegrator.dp_db_gyro());

    dT_ = DualValuef(preintegrator.dT());

    setInformationMatrix(preintegrator.sigma().inverse());
  }

  BOSS_REGISTER_AND_INSTANTIATE(ImuPreintegrationFactorAD)
} // namespace srrg2_solver