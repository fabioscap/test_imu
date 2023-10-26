#include "variable_se3_expmap_right.h"

#include "srrg_solver/solver_core/instance_macros.h"
#include "srrg_solver/solver_core/variable_impl.cpp"
#include <srrg_geometry/geometry3d.h>

namespace srrg2_solver {
  void VariableSE3ExpMapRight::applyPerturbation(const Vector6f& pert) {
    ;
  }

  void VariableSE3ExpMapRightAD::applyPerturbationAD(const ADPerturbationVectorType& pert) {
    const auto& dp = pert.head(3);
    const auto& dr = pert.tail(3);

    auto dR = geometry3d::expMapSO3(static_cast<Vector3_<DualType>>(dr));

    _ad_estimate.linear() *= dR;

    _ad_estimate.translation() += _ad_estimate.linear() * dp;
  }

  INSTANTIATE(VariableSE3ExpMapRight)
  INSTANTIATE(VariableSE3ExpMapRightAD)

} // namespace srrg2_solver