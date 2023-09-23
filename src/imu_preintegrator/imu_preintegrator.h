#include <common/common.h>

namespace test_imu {
  using CovType = core::MatrixN_<float, 15>;

  struct PreintegratedImuMeasurement {
    float dT;

    core::Vector3f delta_p;
    core::Matrix3f delta_R;
    core::Vector3f delta_v;

    CovType sigma;

    // nominal values for the bias
    core::Vector3f bias_acc;
    core::Vector3f bias_gyro;

    // bias correction
    core::Matrix3f dR_db_gyro;
    core::Matrix3f dv_db_acc;
    core::Matrix3f dv_db_gyro;
    core::Matrix3f dp_db_acc;
    core::Matrix3f dp_db_gyro;

    PreintegratedImuMeasurement() = default;
    inline PreintegratedImuMeasurement(const PreintegratedImuMeasurement& other) {
      std::cout << "Sborra culo merda cazzo\n" << std::endl;
    }
  };

  class ImuPreintegrator {
  public:
    using CovType      = Eigen::MatrixXf;
    using CovNoiseType = Eigen::MatrixXf;
    using AType        = Eigen::MatrixXf;
    using BType        = Eigen::MatrixXf;

    void preintegrate(const ImuMeasurement& m, float dt);

    // set the first measurement, the initial covariance and the biases
    void reset();

    const core::Matrix3f& delta_R() const {
      return delta_R_;
    }
    const core::Vector3f& delta_p() const {
      return delta_p_;
    }
    const core::Vector3f& delta_v() const {
      return delta_v_;
    }

    // DEBUG function
    void getPrediction(const core::Isometry3f& Ti,
                       const core::Vector3f& vi,
                       core::Isometry3f& Tf,
                       core::Vector3f& vf) const;

    // produces a PreintegratedImuMeasurement
    // it COPIES
    void getMeasurement(PreintegratedImuMeasurement& preintegrated_measurement) const;

  protected:
    std::vector<ImuMeasurement> measurements_;

    core::Matrix3f delta_R_ = core::Matrix3f::Identity();
    core::Vector3f delta_v_ = core::Vector3f::Zero();
    core::Vector3f delta_p_ = core::Vector3f::Zero();

    // bias estimates at starting step. You need just these for
    // preintegration.
    core::Vector3f bias_acc_  = core::Vector3f::Zero();
    core::Vector3f bias_gyro_ = core::Vector3f::Zero();

    CovType sigma_ = CovType::Zero(15, 15);
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
    // this should remain the same during operation
    CovNoiseType sigma_noise_ = CovNoiseType::Zero(18, 18);

    // we need discretize the covariances
    // see Optimal state estimation 8.1
    // or https://github.com/borglab/gtsam/blob/develop/doc/ImuFactor.pdf
    CovNoiseType scaling_ = CovNoiseType::Zero(18, 18);
    //
  };

} // namespace test_imu