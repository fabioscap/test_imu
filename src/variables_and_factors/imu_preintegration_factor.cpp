#include "imu_preintegration_factor.h"

#include <srrg_solver/solver_core/ad_error_factor_impl.cpp>
#include <srrg_solver/solver_core/error_factor_impl.cpp>

#include <srrg_solver/solver_core/instance_macros.h>

namespace srrg2_solver {
  // TODO can I reuse code

  ImuPreintegrationFactorAD::ADErrorVectorType
  ImuPreintegrationFactorAD::operator()(ImuPreintegrationFactorAD::VariableTupleType& vars) {
    const Isometry3_<DualValuef>& T_i = vars.at<0>()->adEstimate() * offset_;
    const dMatrix3f& R_i              = T_i.linear();
    const dVector3f& p_i              = T_i.translation();
    const dVector3f& v_i              = offset_.linear().transpose() * vars.at<1>()->adEstimate();

    const Isometry3_<DualValuef>& T_j = vars.at<2>()->adEstimate() * offset_;
    const dMatrix3f& R_j              = T_j.linear();
    const dVector3f& p_j              = T_j.translation();
    const dVector3f& v_j              = offset_.linear().transpose() * vars.at<3>()->adEstimate();

    const dVector3f& bias_acc_i  = vars.at<4>()->adEstimate();
    const dVector3f& bias_gyro_i = vars.at<5>()->adEstimate();
    const dVector3f& bias_acc_j  = vars.at<6>()->adEstimate();
    const dVector3f& bias_gyro_j = vars.at<7>()->adEstimate();

    ADErrorVectorType error;

    // bias correction via first order approximation around nominal values
    const dVector3f delta_bacc  = bias_acc_i - bias_acc_nom_;
    const dVector3f delta_bgyro = bias_gyro_i - bias_gyro_nom_;

    // orientation error
    error.segment<3>(0) = geometry3d::logMapSO3(
      ((delta_R_ * geometry3d::expMapSO3((dR_db_gyro_ * delta_bgyro).eval())).transpose() *
       R_i.transpose() * R_j)
        .eval());
    // velocity error
    error.segment<3>(3) = R_i.transpose() * (v_j - v_i - grav_ * dT_) -
                          (delta_v_ + dv_db_acc_ * delta_bacc + dv_db_gyro_ * delta_bgyro);
    // position error
    error.segment<3>(6) =
      R_i.transpose() * (p_j - p_i - v_i * dT_ - DualValuef(0.5) * grav_ * dT_ * dT_) -
      (delta_p_ + dp_db_acc_ * delta_bacc + dp_db_gyro_ * delta_bgyro);

    error.segment<3>(9)  = bias_acc_j - bias_acc_i;
    error.segment<3>(12) = bias_gyro_j - bias_gyro_i;

    /*     std::cout << "error:\n"
                  << error(0).value << " " << error(1).value << " " << error(2).value << " "
                  << error(3).value << " " << error(4).value << " " << error(5).value << " "
                  << error(6).value << " " << error(7).value << " " << error(8).value << " "
                  << "\n";
     */
    return error;
  }

  void ImuPreintegrationFactorAD::setMeasurement(test_imu::ImuPreintegratorBase& preintegrator) {
    convertMatrix(delta_R_, preintegrator.delta_R());
    convertMatrix(delta_v_, preintegrator.delta_v());
    convertMatrix(delta_p_, preintegrator.delta_p());
    convertMatrix(bias_acc_nom_, preintegrator.bias_acc());
    convertMatrix(bias_gyro_nom_, preintegrator.bias_gyro());

    const test_imu::BiasJacobians* bias_J_ptr = preintegrator.biasJacobians();
    if (bias_J_ptr) {
      convertMatrix(dR_db_gyro_, bias_J_ptr->dR_db_gyro);
      convertMatrix(dv_db_acc_, bias_J_ptr->dv_db_acc);
      convertMatrix(dv_db_gyro_, bias_J_ptr->dv_db_gyro);
      convertMatrix(dp_db_acc_, bias_J_ptr->dp_db_acc);
      convertMatrix(dp_db_gyro_, bias_J_ptr->dp_db_gyro);
    }
    dT_ = DualValuef(preintegrator.dT());
    setInformationMatrix(preintegrator.sigma().inverse().cast<float>());
    // setInformationMatrix(Eigen::MatrixXf::Identity(15, 15));
  }

  ImuPreintegrationFactorSlimAD::ADErrorVectorType ImuPreintegrationFactorSlimAD::operator()(
    ImuPreintegrationFactorSlimAD::VariableTupleType& vars) {
    const Isometry3_<DualValuef>& T_i = vars.at<0>()->adEstimate() * offset_;
    const dMatrix3f& R_i              = T_i.linear();
    const dVector3f& p_i              = T_i.translation();
    const dVector3f& v_i              = offset_.linear().transpose() * vars.at<1>()->adEstimate();

    const Isometry3_<DualValuef>& T_j = vars.at<2>()->adEstimate() * offset_;
    const dMatrix3f& R_j              = T_j.linear();
    const dVector3f& p_j              = T_j.translation();
    const dVector3f& v_j              = offset_.linear().transpose() * vars.at<3>()->adEstimate();

    const dVector3f& bias_acc_i  = vars.at<4>()->adEstimate();
    const dVector3f& bias_gyro_i = vars.at<5>()->adEstimate();

    ADErrorVectorType error;

    // bias correction via first order approximation around nominal values
    const dVector3f delta_bacc  = bias_acc_i - bias_acc_nom_;
    const dVector3f delta_bgyro = bias_gyro_i - bias_gyro_nom_;

    // orientation error
    error.segment<3>(0) = geometry3d::logMapSO3(
      ((delta_R_ * geometry3d::expMapSO3((dR_db_gyro_ * delta_bgyro).eval())).transpose() *
       R_i.transpose() * R_j)
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

  void
  ImuPreintegrationFactorSlimAD::setMeasurement(test_imu::ImuPreintegratorBase& preintegrator) {
    convertMatrix(delta_R_, preintegrator.delta_R());
    convertMatrix(delta_v_, preintegrator.delta_v());
    convertMatrix(delta_p_, preintegrator.delta_p());
    convertMatrix(bias_acc_nom_, preintegrator.bias_acc());
    convertMatrix(bias_gyro_nom_, preintegrator.bias_gyro());

    const test_imu::BiasJacobians* bias_J_ptr = preintegrator.biasJacobians();
    if (bias_J_ptr) {
      convertMatrix(dR_db_gyro_, bias_J_ptr->dR_db_gyro);
      convertMatrix(dv_db_acc_, bias_J_ptr->dv_db_acc);
      convertMatrix(dv_db_gyro_, bias_J_ptr->dv_db_gyro);
      convertMatrix(dp_db_acc_, bias_J_ptr->dp_db_acc);
      convertMatrix(dp_db_gyro_, bias_J_ptr->dp_db_gyro);
    }
    dT_ = DualValuef(preintegrator.dT());
    setInformationMatrix(preintegrator.sigma().inverse().cast<float>());
    // setInformationMatrix(Eigen::MatrixXf::Identity(9, 9));
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

  void ImuPreintegrationFactor::errorAndJacobian(bool error_only_) {
    const Isometry3f& T_i = _variables.at<0>()->estimate() * offset_;
    const Matrix3f& R_i   = T_i.linear();
    const Vector3f& p_i   = T_i.translation();
    const Vector3f& v_i   = offset_.linear().transpose() * _variables.at<1>()->estimate();

    const Isometry3f& T_j = _variables.at<2>()->estimate() * offset_;
    const Matrix3f& R_j   = T_j.linear();
    const Vector3f& p_j   = T_j.translation();
    const Vector3f& v_j   = offset_.linear().transpose() * _variables.at<3>()->estimate();

    const Vector3f& bias_acc_i  = _variables.at<4>()->estimate();
    const Vector3f& bias_gyro_i = _variables.at<5>()->estimate();
    const Vector3f& bias_acc_j  = _variables.at<6>()->estimate();
    const Vector3f& bias_gyro_j = _variables.at<7>()->estimate();

    // bias correction via first order approximation around nominal values
    const Vector3f delta_bacc  = bias_acc_i - bias_acc_nom_;
    const Vector3f delta_bgyro = bias_gyro_i - bias_gyro_nom_;

    _e.segment<3>(0) = geometry3d::logMapSO3(
      ((delta_R_ * geometry3d::expMapSO3((dR_db_gyro_ * delta_bgyro).eval())).transpose() *
       R_i.transpose() * R_j)
        .eval());
    // velocity error
    _e.segment<3>(3) = R_i.transpose() * (v_j - v_i - grav_ * dT_) -
                       (delta_v_ + dv_db_acc_ * delta_bacc + dv_db_gyro_ * delta_bgyro);
    // position error
    _e.segment<3>(6) =
      R_i.transpose() * (p_j - p_i - v_i * dT_ - DualValuef(0.5) * grav_ * dT_ * dT_) -
      (delta_p_ + dp_db_acc_ * delta_bacc + dp_db_gyro_ * delta_bgyro);

    _e.segment<3>(9)  = (bias_acc_j - bias_acc_i);
    _e.segment<3>(12) = (bias_gyro_j - bias_gyro_i);

    if (error_only_)
      return;

    // std::cout << "error: " << _e.transpose() << "\n";

    // check if this is necessary
    // std::cout <<jacobian<0>() << "\n";
    jacobian<0>().setZero();
    jacobian<1>().setZero();
    jacobian<2>().setZero();
    jacobian<3>().setZero();
    jacobian<4>().setZero();
    jacobian<5>().setZero();
    jacobian<6>().setZero();
    jacobian<7>().setZero();

    // jacobians of delta_p_
    auto Jpi  = jacobian<0>(); // 15x6
    auto Jvi  = jacobian<1>(); // 15x3
    auto Jpj  = jacobian<2>(); // 15x6
    auto Jvj  = jacobian<3>(); // 15x3
    auto Jbai = jacobian<4>(); // 15x3
    auto Jbgi = jacobian<5>(); // 15x3
    auto Jbaj = jacobian<6>(); // 15x3
    auto Jbgj = jacobian<7>(); // 15x3

    // jacobians related to delta_p_, which starts at the 6th row
    Jpi.block<3, 3>(6, 3) = geometry3d::skew(
      (R_i.transpose() * (p_j - p_i - v_i * dT_ - 0.5 * grav_ * dT_ * dT_)).eval());

    Jpi.block<3, 3>(6, 0) = -Matrix3f::Identity();

    Jvi.block<3, 3>(6, 0) = -R_i.transpose() * dT_;

    Jbai.block<3, 3>(6, 0) = -dp_db_acc_;

    Jpj.block<3, 3>(6, 0) = R_i.transpose() * R_j;

    Jbgi.block<3, 3>(6, 0) = -dp_db_gyro_;

    // jacobians related to delta_v_
    Jpi.block<3, 3>(3, 3) = geometry3d::skew((R_i.transpose() * (v_j - v_i - grav_ * dT_)).eval());

    Jvi.block<3, 3>(3, 0) = -R_i.transpose();

    Jbai.block<3, 3>(3, 0) = -dv_db_acc_;

    Jvj.block<3, 3>(3, 0) = R_i.transpose();

    Jbgi.block<3, 3>(3, 0) = -dv_db_gyro_;

    // jacobians related to delta_R_
    Vector3f r_phi        = _e.segment<3>(0);
    Matrix3f Jinv         = jacobianExpMapSO3inv(r_phi);
    Jpi.block<3, 3>(0, 3) = -Jinv * R_j.transpose() * R_i;

    Jpj.block<3, 3>(0, 3) = Jinv;

    Jbgi.block<3, 3>(0, 0) = -Jinv * geometry3d::expMapSO3(r_phi).transpose() *
                             geometry3d::jacobianExpMapSO3((dR_db_gyro_ * delta_bgyro).eval()) *
                             dR_db_gyro_;

    // now the jacobians related to the bias propagation
    Jbaj.block<3, 3>(9, 0)  = Matrix3f::Identity();
    Jbgj.block<3, 3>(12, 0) = Matrix3f::Identity();

    Jbai.block<3, 3>(9, 0)  = -Matrix3f::Identity();
    Jbgi.block<3, 3>(12, 0) = -Matrix3f::Identity();

    // std::cout << _J << "\n\n";

    /* std::cout << "error: " << _e.transpose() << "\n";
    std::cout << "J:\n" << _J << "\n\n\n"; */
    return;
  }

  void ImuPreintegrationFactorSlim::errorAndJacobian(bool error_only_) {
    const Isometry3f& T_i = _variables.at<0>()->estimate() * offset_;
    const Matrix3f& R_i   = T_i.linear();
    const Vector3f& p_i   = T_i.translation();
    const Vector3f& v_i   = offset_.linear().transpose() * _variables.at<1>()->estimate();

    const Isometry3f& T_j = _variables.at<2>()->estimate() * offset_;
    const Matrix3f& R_j   = T_j.linear();
    const Vector3f& p_j   = T_j.translation();
    const Vector3f& v_j   = offset_.linear().transpose() * _variables.at<3>()->estimate();

    const Vector3f& bias_acc_i  = _variables.at<4>()->estimate();
    const Vector3f& bias_gyro_i = _variables.at<5>()->estimate();

    // bias correction via first order approximation around nominal values
    const Vector3f delta_bacc  = bias_acc_i - bias_acc_nom_;
    const Vector3f delta_bgyro = bias_gyro_i - bias_gyro_nom_;

    _e.segment<3>(0) = geometry3d::logMapSO3(
      ((delta_R_ * geometry3d::expMapSO3((dR_db_gyro_ * delta_bgyro).eval())).transpose() *
       R_i.transpose() * R_j)
        .eval());
    // velocity error
    _e.segment<3>(3) = R_i.transpose() * (v_j - v_i - grav_ * dT_) -
                       (delta_v_ + dv_db_acc_ * delta_bacc + dv_db_gyro_ * delta_bgyro);
    // position error
    _e.segment<3>(6) =
      R_i.transpose() * (p_j - p_i - v_i * dT_ - DualValuef(0.5) * grav_ * dT_ * dT_) -
      (delta_p_ + dp_db_acc_ * delta_bacc + dp_db_gyro_ * delta_bgyro);

    if (error_only_)
      return;

    // std::cout << "error: " << _e.transpose() << "\n";

    // check if this is necessary
    // std::cout <<jacobian<0>() << "\n";
    jacobian<0>().setZero();
    jacobian<1>().setZero();
    jacobian<2>().setZero();
    jacobian<3>().setZero();
    jacobian<4>().setZero();
    jacobian<5>().setZero();

    // jacobians of delta_p_
    auto Jpi  = jacobian<0>(); // 15x6
    auto Jvi  = jacobian<1>(); // 15x3
    auto Jpj  = jacobian<2>(); // 15x6
    auto Jvj  = jacobian<3>(); // 15x3
    auto Jbai = jacobian<4>(); // 15x3
    auto Jbgi = jacobian<5>(); // 15x3

    // jacobians related to delta_p_, which starts at the 6th row
    Jpi.block<3, 3>(6, 3) = geometry3d::skew(
      (R_i.transpose() * (p_j - p_i - v_i * dT_ - 0.5 * grav_ * dT_ * dT_)).eval());

    Jpi.block<3, 3>(6, 0) = -Matrix3f::Identity();

    Jvi.block<3, 3>(6, 0) = -R_i.transpose() * dT_;

    Jbai.block<3, 3>(6, 0) = -dp_db_acc_;

    Jpj.block<3, 3>(6, 0) = R_i.transpose() * R_j;

    Jbgi.block<3, 3>(6, 0) = -dp_db_gyro_;

    // jacobians related to delta_v_
    Jpi.block<3, 3>(3, 3) = geometry3d::skew((R_i.transpose() * (v_j - v_i - grav_ * dT_)).eval());

    Jvi.block<3, 3>(3, 0) = -R_i.transpose();

    Jbai.block<3, 3>(3, 0) = -dv_db_acc_;

    Jvj.block<3, 3>(3, 0) = R_i.transpose();

    Jbgi.block<3, 3>(3, 0) = -dv_db_gyro_;

    // jacobians related to delta_R_
    Vector3f r_phi        = _e.segment<3>(0);
    Matrix3f Jinv         = jacobianExpMapSO3inv(r_phi);
    Jpi.block<3, 3>(0, 3) = -Jinv * R_j.transpose() * R_i;

    Jpj.block<3, 3>(0, 3) = Jinv;

    Jbgi.block<3, 3>(0, 0) = -Jinv * geometry3d::expMapSO3(r_phi).transpose() *
                             geometry3d::jacobianExpMapSO3((dR_db_gyro_ * delta_bgyro).eval()) *
                             dR_db_gyro_;

    /* std::cout << "error: " << _e.transpose() << "\n";
    std::cout << "J:\n" << _J << "\n\n\n"; */

    return;
  }

  void ImuPreintegrationFactor::setMeasurement(test_imu::ImuPreintegratorBase& preintegrator) {
    (delta_R_ = preintegrator.delta_R());
    (delta_v_ = preintegrator.delta_v());
    (delta_p_ = preintegrator.delta_p());
    (bias_acc_nom_ = preintegrator.bias_acc());
    (bias_gyro_nom_ = preintegrator.bias_gyro());

    const test_imu::BiasJacobians* bias_J_ptr = preintegrator.biasJacobians();
    if (bias_J_ptr) {
      (dR_db_gyro_ = bias_J_ptr->dR_db_gyro.cast<float>());
      (dv_db_acc_ = bias_J_ptr->dv_db_acc.cast<float>());
      (dv_db_gyro_ = bias_J_ptr->dv_db_gyro.cast<float>());
      (dp_db_acc_ = bias_J_ptr->dp_db_acc.cast<float>());
      (dp_db_gyro_ = bias_J_ptr->dp_db_gyro.cast<float>());
    }
    dT_ = preintegrator.dT();
    setInformationMatrix(preintegrator.sigma().inverse().cast<float>());
    // setInformationMatrix(Eigen::MatrixXf::Identity(15, 15));
  }

  void ImuPreintegrationFactorSlim::setMeasurement(test_imu::ImuPreintegratorBase& preintegrator) {
    (delta_R_ = preintegrator.delta_R());
    (delta_v_ = preintegrator.delta_v());
    (delta_p_ = preintegrator.delta_p());
    (bias_acc_nom_ = preintegrator.bias_acc());
    (bias_gyro_nom_ = preintegrator.bias_gyro());

    const test_imu::BiasJacobians* bias_J_ptr = preintegrator.biasJacobians();
    if (bias_J_ptr) {
      (dR_db_gyro_ = bias_J_ptr->dR_db_gyro.cast<float>());
      (dv_db_acc_ = bias_J_ptr->dv_db_acc.cast<float>());
      (dv_db_gyro_ = bias_J_ptr->dv_db_gyro.cast<float>());
      (dp_db_acc_ = bias_J_ptr->dp_db_acc.cast<float>());
      (dp_db_gyro_ = bias_J_ptr->dp_db_gyro.cast<float>());
    }
    dT_ = preintegrator.dT();
    setInformationMatrix(preintegrator.sigma().inverse().cast<float>());
    // setInformationMatrix(Eigen::MatrixXf::Identity(9, 9));
  }

  void ImuPreintegrationFactor::_drawImpl(ViewerCanvasPtr canvas_) const {
    Vector3f coords[2];
    coords[0] =
      reinterpret_cast<const VariableSE3ExpMapRight*>(variable(0))->estimate().translation();
    coords[1] =
      reinterpret_cast<const VariableSE3ExpMapRight*>(variable(2))->estimate().translation();

    float lw = 0.5;
    if (fabs(variableId(0) - variableId(1)) == 1) {
      lw *= 2;
    }
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

  void ImuPreintegrationFactorAD::_drawImpl(ViewerCanvasPtr canvas_) const {
    Vector3f coords[2];
    coords[0] =
      reinterpret_cast<const VariableSE3ExpMapRight*>(variable(0))->estimate().translation();
    coords[1] =
      reinterpret_cast<const VariableSE3ExpMapRight*>(variable(2))->estimate().translation();

    float lw = 0.5;
    if (fabs(variableId(0) - variableId(1)) == 1) {
      lw *= 2;
    }
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

  void ImuPreintegrationFactorSlim::_drawImpl(ViewerCanvasPtr canvas_) const {
    Vector3f coords[2];
    coords[0] =
      reinterpret_cast<const VariableSE3ExpMapRight*>(variable(0))->estimate().translation();
    coords[1] =
      reinterpret_cast<const VariableSE3ExpMapRight*>(variable(2))->estimate().translation();

    float lw = 0.5;
    if (fabs(variableId(0) - variableId(1)) == 1) {
      lw *= 2;
    }
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

  void ImuPreintegrationFactorSlimAD::_drawImpl(ViewerCanvasPtr canvas_) const {
    Vector3f coords[2];
    coords[0] =
      reinterpret_cast<const VariableSE3ExpMapRight*>(variable(0))->estimate().translation();
    coords[1] =
      reinterpret_cast<const VariableSE3ExpMapRight*>(variable(2))->estimate().translation();

    float lw = 0.5;
    if (fabs(variableId(0) - variableId(1)) == 1) {
      lw *= 2;
    }
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

  INSTANTIATE(ImuPreintegrationFactor)
  INSTANTIATE(ImuPreintegrationFactorAD)
  INSTANTIATE(ImuPreintegrationFactorSlim)
  INSTANTIATE(ImuPreintegrationFactorSlimAD) INSTANTIATE(BiasErrorFactorAD)

} // namespace srrg2_solver