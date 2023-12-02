#include "gps_factor_ad.h"
#include <srrg_solver/solver_core/ad_error_factor_impl.cpp>
#include <srrg_solver/solver_core/error_factor_impl.cpp>
#include <srrg_solver/solver_core/instance_macros.h>

namespace srrg2_solver {

  void GpsFactorAD::_drawImpl(ViewerCanvasPtr canvas_) const {
    if (!canvas_) {
      throw std::runtime_error("VariablePoint3::draw|invalid canvas");
    }
    canvas_->pushColor();
    canvas_->setColor(srrg2_core::ColorPalette::color3fCyan());
    canvas_->pushMatrix();
    Eigen::Isometry3f iso = Eigen::Isometry3f::Identity();
    iso.translation()     = _measurement;
    canvas_->multMatrix(iso.matrix());
    canvas_->putSphere(0.1);
    canvas_->popMatrix();
    canvas_->popAttribute();

    Vector3f coords[2];
    coords[0] = _measurement;
    coords[1] =
      reinterpret_cast<const VariableSE3ExpMapRight*>(variable(0))->estimate().translation();

    float lw = 0.5;
    lw *= (level() * 3 + 1);
    canvas_->pushColor();
    canvas_->pushLineWidth();
    canvas_->setLineWidth(lw);
    float fading   = 1. - 0.5 * level();
    Vector3f color = srrg2_core::ColorPalette::color3fBlue() * fading;
    canvas_->setColor(color);
    canvas_->putLine(2, coords);
    canvas_->popAttribute();
    canvas_->popAttribute();
  }

  INSTANTIATE(GpsFactorAD)
} // namespace srrg2_solver
