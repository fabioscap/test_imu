//! include all the types you declared here
#include "gps_factor_ad.h"
#include "imu_preintegration_factor.h"
#include "variable_se3_expmap_right.h"

//! we try to instantiate all solvers
namespace srrg2_solver {
  using namespace srrg2_core;

  // this is the function you have to call to initialize
  // the serialization subsystem
  void variables_and_factors_imu_registerTypes() {
    BOSS_REGISTER_CLASS(VariableSE3ExpMapRight)
    BOSS_REGISTER_CLASS(VariableSE3ExpMapRightAD)
    BOSS_REGISTER_CLASS(ImuPreintegrationFactorAD)
    BOSS_REGISTER_CLASS(ImuPreintegrationFactorSlimAD)
    BOSS_REGISTER_CLASS(BiasErrorFactorAD)
    BOSS_REGISTER_CLASS(GpsFactorAD)
  }
} // namespace srrg2_solver
