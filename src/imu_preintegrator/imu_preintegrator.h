#include <common/common.h>

namespace test_imu {
  class ImuPreintegrator {
  public:
    using CovType      = core::MatrixN_<float, 15>;
    using CovNoiseType = Eigen::MatrixXf;
    using AType        = Eigen::MatrixXf;
    using BType        = Eigen::MatrixXf;
    // perform preintegration up until m_new.timestamp
    // therefore the measurement m_new is used in the next
    // call to preintegrate
    // this is a bit backward
    // to change this add a parameter dt integration step and perform operations on m_new
    void preintegrate(const ImuMeasurement& m_new);

    // set the first measurement, the initial covariance and the biases
    void reset(const ImuMeasurement& measurement);

    const core::Matrix3f& delta_R() const {
      return delta_R_;
    }
    const core::Vector3f& delta_p() const {
      return delta_p_;
    }
    const core::Vector3f& delta_v() const {
      return delta_v_;
    }

    void getPrediction(const core::Isometry3f& Ti,
                       const core::Vector3f& vi,
                       core::Isometry3f& Tf,
                       core::Vector3f& vf) const;

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
    // Eigen::Block<CovType, 3, 3> sigma_R_         = sigma_.block<3, 3>(3, 3);
    // Eigen::Block<CovType, 3, 3> sigma_p_         = sigma_.block<3, 3>(3, 3);
    // Eigen::Block<CovType, 3, 3> sigma_v_         = sigma_.block<3, 3>(3, 3);
    // Eigen::Block<CovType, 3, 3> sigma_bias_acc_  = sigma_.block<3, 3>(3, 3);
    // Eigen::Block<CovType, 3, 3> sigma_bias_gryo_ = sigma_.block<3, 3>(3, 3);

    // sensor-dependent
    core::Matrix6f sigma_imu_ = core::Matrix6f::Zero();

    float t_ = 0; // preintegration time

    // store matrices for noise propagation
    AType A_ = AType::Identity(15, 15);
    BType B_ = BType::Zero(15, 18);

    // store matrices for bias correction
    core::Matrix3f dR_db_gyro_ = core::Matrix3f::Zero();
    core::Matrix3f dv_db_acc_  = core::Matrix3f::Zero();
    core::Matrix3f dv_db_gyro_ = core::Matrix3f::Zero();
    core::Matrix3f dp_db_acc_  = core::Matrix3f::Zero();
    core::Matrix3f dp_db_gyro_ = core::Matrix3f::Zero();

    // block diagonal
    CovNoiseType sigma_noise_ = CovNoiseType::Zero(18, 18);
    //
  };

} // namespace test_imu