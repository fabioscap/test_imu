#include <common/common.h>

#include <srrg_solver/solver_core/error_factor.h>
#include <srrg_solver/variables_and_factors/types_3d/variable_se3.h>
#include <srrg_solver/variables_and_factors/types_common/all_types.h>

namespace test_imu {

  namespace solver = srrg2_solver;

  class ImuPreintegrationFactor
    : public solver::ErrorFactor_<15,
                                  solver::VariableSE3QuaternionRight, // pose_from
                                  solver::VariableVector3,            // vel_from
                                  solver::VariableSE3QuaternionRight, // pose_to
                                  solver::VariableVector3,            // vel_to
                                  solver::VariableVector3,            // bias_acc
                                  solver::VariableVector3> {          // bias_gyro
  public:
    // get variable: _variables.at<idx>()
    // store error in _e
    // store jacobian de/dvidx in jacobian<idx>()

    // see if I can subclass pose-pose constraint

    void errorAndJacobian(bool error_only) override;

  protected:
  };

  class ImuPreintegrator {
  public:
    using CovType = core::MatrixN_<float, 15>;
    using AType   = core::MatrixN_<float, 9>;
    using BType   = Eigen::Matrix<float, 9, 6>;

    void preintegrate(const ImuMeasurement& measurement);

    // set the first measurement, the initial covariance and the biases
    void reset(const ImuMeasurement& measurement);

  protected:
    std::vector<ImuMeasurement> measurements_;

    core::Matrix3f delta_R_ = core::Matrix3f::Identity();
    core::Vector3f delta_v_ = core::Vector3f::Zero();
    core::Vector3f delta_p_ = core::Vector3f::Zero();

    // bias estimates at starting step. You need just these for
    // preintegration.
    core::Vector3f bias_acc_  = core::Vector3f::Zero();
    core::Vector3f bias_gyro_ = core::Vector3f::Zero();

    CovType sigma_;
    Eigen::Block<CovType, 3, 3> sigma_R_         = sigma_.block<3, 3>(3, 3);
    Eigen::Block<CovType, 3, 3> sigma_p_         = sigma_.block<3, 3>(3, 3);
    Eigen::Block<CovType, 3, 3> sigma_v_         = sigma_.block<3, 3>(3, 3);
    Eigen::Block<CovType, 3, 3> sigma_bias_acc_  = sigma_.block<3, 3>(3, 3);
    Eigen::Block<CovType, 3, 3> sigma_bias_gryo_ = sigma_.block<3, 3>(3, 3);

    // sensor-dependent
    core::Matrix6f sigma_imu_ = core::Matrix6f::Zero();

    float t_ = 0; // preintegration time

    // store matrices for noise propagation
    AType A_ = AType::Identity();
    BType B_ = BType::Zero();
  };

} // namespace test_imu