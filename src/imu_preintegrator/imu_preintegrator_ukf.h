#include "common/common.h"
#include "common/manifold.h"
#include "common/unscented.h"

namespace test_imu {

  class ImuPreintegratorUKF {
  public:
    using CovType      = Eigen::MatrixXf;
    using CovNoiseType = Eigen::MatrixXf;

    static constexpr int state_dim = 15;
    static constexpr int noise_dim = 12;

    using CovJType = core::MatrixN_<float, state_dim + noise_dim>;

    void preintegrate(const ImuMeasurement& m, float dt);
    // allocates a new imu measurement
    void reset();
    // clang-format off
    inline const core::Matrix3f& delta_R() const {return delta_incr_.get<0>().data();}
    inline const core::Vector3f& delta_p() const {return delta_incr_.get<2>().data();}
    inline const core::Vector3f& delta_v() const {return delta_incr_.get<1>().data();}
    inline const CovType& sigma() {
      sigma_ = sigma_joint_.block<state_dim, state_dim>(0, 0);
      return sigma_;}
    inline const core::Vector3f& bias_acc() const {return bias_acc_;}
    inline const core::Vector3f& bias_gyro() const {return bias_gyro_;}
    inline const core::Matrix3f& dR_db_gyro() const { return dR_db_gyro_;}
    inline const core::Matrix3f& dv_db_acc() const { return dv_db_acc_;}
    inline const core::Matrix3f& dv_db_gyro() const { return dv_db_gyro_;}
    inline const core::Matrix3f& dp_db_acc() const { return dp_db_acc_;}
    inline const core::Matrix3f& dp_db_gyro() const { return dp_db_gyro_;}
    inline const float dT() const {return dT_;}
    // clang-format on

    const void setNoiseGyroscope(const core::Vector3f& v);
    const void setNoiseAccelerometer(const core::Vector3f& v);
    const void setNoiseBiasGyroscope(const core::Vector3f& v);
    const void setNoiseBiasAccelerometer(const core::Vector3f& v);

    using DeltaManifold = ManifoldComp_<ManifoldSO3,   // dR
                                        Euclidean_<3>, // dv
                                        Euclidean_<3>, // dp
                                        Euclidean_<3>, // ba
                                        Euclidean_<3>, // bg
                                        Euclidean_<3>, // na
                                        Euclidean_<3>, // ng
                                        Euclidean_<3>, // nba
                                        Euclidean_<3>  // nbg
                                        >;
    // DEBUG function
    void getPrediction(const core::Isometry3f& Ti,
                       const core::Vector3f& vi,
                       core::Isometry3f& Tf,
                       core::Vector3f& vf);

    // protected:
    std::vector<ImuMeasurement> measurements_;

    DeltaManifold delta_incr_;

    float dT_ = 0;

    // UKF: preallocate sigma joint input and state
    CovJType sigma_joint_ = 1e-6 * CovJType::Identity();
    CovType sigma_;

    // container for sigma points
    SigmaPoints<DeltaManifold> spoints;

    // nominal values for the bias
    core::Vector3f bias_acc_  = core::Vector3f::Zero();
    core::Vector3f bias_gyro_ = core::Vector3f::Zero();

    // bias correction
    core::Matrix3f dR_db_gyro_ = core::Matrix3f::Zero();
    core::Matrix3f dv_db_acc_  = core::Matrix3f::Zero();
    core::Matrix3f dv_db_gyro_ = core::Matrix3f::Zero();
    core::Matrix3f dp_db_acc_  = core::Matrix3f::Zero();
    core::Matrix3f dp_db_gyro_ = core::Matrix3f::Zero();

    // this should remain the same during operation
    CovNoiseType sigma_noise_ = 1e-9 * CovNoiseType::Identity(noise_dim, noise_dim);

    // we need discretize the covariances
    // see Optimal state estimation 8.1
    // or https://github.com/borglab/gtsam/blob/develop/doc/ImuFactor.pdf
    CovNoiseType scaling_ = CovNoiseType::Identity(noise_dim, noise_dim);
    //
  };

} // namespace test_imu