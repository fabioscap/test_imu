#include <common/common.h>

namespace test_imu {

  // this works but the covariance is always singular, even after just 1 step
  // I believe this is due to the fact that the process becomes unobservable when adding also the
  // biases but it could also be bugs in covariance propagation
  // solved by removing initial uncertainty on biases
  class ImuPreintegrator {
  public:
    using CovType      = Eigen::MatrixXf;
    using CovNoiseType = Eigen::MatrixXf;
    using AType        = Eigen::MatrixXf;
    using BType        = Eigen::MatrixXf;

    static constexpr int state_dim = 15;
    static constexpr int noise_dim = 12;

    void preintegrate(const ImuMeasurement& m, float dt);
    // allocates a new imu measurement
    void reset();
    // clang-format off
    inline const core::Matrix3f& delta_R() const {return delta_R_;}
    inline const core::Vector3f& delta_p() const {return delta_p_;}
    inline const core::Vector3f& delta_v() const {return delta_v_;}
    inline const CovType& sigma() const {return sigma_;}
    inline const core::Vector3f& bias_acc() const {return bias_acc_;}
    inline const core::Vector3f& bias_gyro() const {return bias_gyro_;}
    inline const core::Matrix3f& dR_db_gyro() const { return dR_db_gyro_;}
    inline const core::Matrix3f& dv_db_acc() const { return dv_db_acc_;}
    inline const core::Matrix3f& dv_db_gyro() const { return dv_db_gyro_;}
    inline const core::Matrix3f& dp_db_acc() const { return dp_db_acc_;}
    inline const core::Matrix3f& dp_db_gyro() const { return dp_db_gyro_;}
    inline const float dT() const {return dT_;}
    // clang-format on

    // DEBUG function
    void getPrediction(const core::Isometry3f& Ti,
                       const core::Vector3f& vi,
                       core::Isometry3f& Tf,
                       core::Vector3f& vf);

    // protected:
    std::vector<ImuMeasurement> measurements_;

    float dT_ = 0;

    core::Vector3f delta_p_ = core::Vector3f::Zero();
    core::Matrix3f delta_R_ = core::Matrix3f::Identity();
    core::Vector3f delta_v_ = core::Vector3f::Zero();

    CovType sigma_ = CovType::Zero(state_dim, state_dim);

    // nominal values for the bias
    core::Vector3f bias_acc_  = core::Vector3f::Zero();
    core::Vector3f bias_gyro_ = core::Vector3f::Zero();

    // bias correction
    core::Matrix3f dR_db_gyro_ = core::Matrix3f::Zero();
    core::Matrix3f dv_db_acc_  = core::Matrix3f::Zero();
    core::Matrix3f dv_db_gyro_ = core::Matrix3f::Zero();
    core::Matrix3f dp_db_acc_  = core::Matrix3f::Zero();
    core::Matrix3f dp_db_gyro_ = core::Matrix3f::Zero();

    // allocate matrices for noise propagation
    AType A_ = AType::Identity(state_dim, state_dim);
    BType B_ = BType::Zero(state_dim, noise_dim);

    // block diagonal
    // this should remain the same during operation
    CovNoiseType sigma_noise_ = 1e-3 * CovNoiseType::Identity(noise_dim, noise_dim);

    // we need discretize the covariances
    // see Optimal state estimation 8.1
    // or https://github.com/borglab/gtsam/blob/develop/doc/ImuFactor.pdf
    CovNoiseType scaling_ = CovNoiseType::Identity(noise_dim, noise_dim);
    //
  };

} // namespace test_imu