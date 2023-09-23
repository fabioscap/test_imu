#include "imu_preintegration_factor.h"

// #include <srrg_solver/solver_core/ad_error_factor_impl.cpp>
// #include <srrg_solver/solver_core/error_factor_impl.cpp>
#include <srrg_solver/solver_core/instance_macros.h>

namespace srrg2_solver {
  ImuPreintegrationFactorAD::ADErrorVectorType
  ImuPreintegrationFactorAD::operator()(ImuPreintegrationFactorAD::VariableTupleType& vars) {
    const Isometry3_<DualValuef>& Ti = vars.at<0>()->adEstimate();
    const Vector3_<DualValuef>& vi   = vars.at<1>()->adEstimate();

    const Isometry3_<DualValuef>& Tj = vars.at<2>()->adEstimate();
    const Vector3_<DualValuef>& vj   = vars.at<3>()->adEstimate();

    const Vector3_<DualValuef>& bias_acc_i  = vars.at<4>()->adEstimate();
    const Vector3_<DualValuef>& bias_acc_j  = vars.at<5>()->adEstimate();
    const Vector3_<DualValuef>& bias_gyro_i = vars.at<6>()->adEstimate();
    const Vector3_<DualValuef>& bias_gyro_j = vars.at<7>()->adEstimate();
    return ImuPreintegrationFactorAD::ADErrorVectorType();
  }

  // INSTANTIATE(ImuPreintegrationFactorAD)
  BOSS_REGISTER_AND_INSTANTIATE(ImuPreintegrationFactorAD)
} // namespace srrg2_solver