#include "imu_preintegration_factor.h"

// #include <srrg_solver/solver_core/ad_error_factor_impl.cpp>
// #include <srrg_solver/solver_core/error_factor_impl.cpp>
#include <srrg_solver/solver_core/instance_macros.h>

namespace srrg2_solver {
  ImuPreintegrationFactorAD::ADErrorVectorType
  ImuPreintegrationFactorAD::operator()(ImuPreintegrationFactorAD::VariableTupleType& vars) {
    using dMatrix3f Matrix3_<DualValuef>;
    using dVector3f Vector3_<DualValuef>;

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

    const dVector3f delta_bacc  = bias_acc_i - pm_.bias_acc;
    const dVector3f delta_bgyro = bias_gyro_i - pm_bias_gyro;

    // we can also use any pose-pose error in place of the one in the paper

    // orientation error
    error.segment<3>(0) = geometry3d::logMapSO3(
      pm_.deltaR * geometry3d::expMapSO3(pm_.dR_db_gyro * pm_.delta_bgyro) * R_i.transpose() * R_j);
    // velocity error
    error.segment<3>(3) =
      R_i.transpose() * (v_j - v_i - grav_ * pm_.dT) -
      (pm_.delta_v + pm_.d_v_db_acc * delta_bacc + pm_.d_v_db_gyro * delta_bgyro);
    // position error
    error.segment<3>(6) =
      R_i.transpose() * (p_j - p_i - v_i * pm_.dT - 0.5 * grav_ * pm_.dT * pm_.dT) -
      (pm_.delta_p + pm_.d_p_db_acc * delta_bacc + pm_.d_p_db_gyro * delta_bgyro);

    return ADErrorVectorType();
  }

  // INSTANTIATE(ImuPreintegrationFactorAD)
  BOSS_REGISTER_AND_INSTANTIATE(ImuPreintegrationFactorAD)
} // namespace srrg2_solver