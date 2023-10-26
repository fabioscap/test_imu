#pragma once

#include <srrg_solver/solver_core/ad_variable.h>

#include "srrg_solver/solver_core/ad_variable.h"
#include <srrg_solver/variables_and_factors/types_3d/variable_se3.h>
// #include <srrg_solver/variables_and_factors/types_3d/variable_se3_ad.h>

namespace srrg2_solver {
  // we implement the same perturbation as in original implementation
  // to have a coherent covariance on chart
  class VariableSE3ExpMapRight : public VariableSE3Base {
  public:
    void applyPerturbation(const Vector6f& pert) override;
  };

  class VariableSE3ExpMapRightAD : public ADVariable_<VariableSE3ExpMapRight> {
  public:
    using DualType                 = ad::DualValue_<float>;
    using VariableType             = VariableSE3ExpMapRight;
    using ADVariableType           = ADVariable_<VariableType>;
    using ADPerturbationVectorType = typename ADVariableType::ADPerturbationVectorType;
    using ADEstimateType           = typename ADVariableType::ADEstimateType;

    void applyPerturbationAD(const ADPerturbationVectorType& pert);
  };
} // namespace srrg2_solver