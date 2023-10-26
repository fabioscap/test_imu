#include <srrg_geometry/geometry_defs.h>

#include "common/common.h"

#include <ctime>
#include <memory>
#include <random>

// generate an SE(3) trajectory
// and produce imu measurements

namespace test_imu {
  // the trajectories must have positive velocity always
  class SE3PlanarTrajectory {
  public:
    SE3PlanarTrajectory(float T) : T_(T){};

    virtual void sampleTrajectory(float t,
                                  core::Vector3f& pos,
                                  core::Vector3f& vel,
                                  core::Vector3f& acc) const {};

    // The frame is always tangent to the curve with its x-axis,
    // and the z-axis always points up (you can add rotation if needed).
    void getPoseMeasurement(float t,
                            core::Isometry3f& pose,
                            core::Vector3f& vel,
                            ImuMeasurement& measurement);

    inline void setT(float T) {
      T_ = T;
    }

    inline float T() const {
      return T_;
    }

  protected:
    float T_; // The duration of the trajectory

    // trajectory transformation parameters
    core::Matrix3f R_ = Rx(0.0f);
    core::Vector3f t_ = core::Vector3f(0, 0, 0);

    // body frame transformation parameters
    core::Matrix3f R_b_ = Ry(0.0f);
  };

  class SE3EightTrajectory : public SE3PlanarTrajectory {
  public:
    SE3EightTrajectory(float T) : SE3PlanarTrajectory(T) {
    }

    void sampleTrajectory(float t,
                          core::Vector3f& pos,
                          core::Vector3f& vel,
                          core::Vector3f& acc) const override;
  };
  class SE3CircleTrajectory : public SE3PlanarTrajectory {
  public:
    SE3CircleTrajectory(float T) : SE3PlanarTrajectory(T) {
    }

    void sampleTrajectory(float t,
                          core::Vector3f& pos,
                          core::Vector3f& vel,
                          core::Vector3f& acc) const override;
  };

  class SE3StraightTrajectory : public SE3PlanarTrajectory {
  public:
    SE3StraightTrajectory(float T) : SE3PlanarTrajectory(T) {
    }

    void sampleTrajectory(float t,
                          core::Vector3f& pos,
                          core::Vector3f& vel,
                          core::Vector3f& acc) const override;
  };

  class FakeImu {
  public:
    FakeImu(std::shared_ptr<SE3PlanarTrajectory> traj_, float freq, int seed = std::time(0)) :
      freq_(freq), trajectory_(traj_), rnd_gen_(seed) {
      core::Vector3f pos, vel, acc;
    }

    void
    generateData(std::vector<std::tuple<ImuMeasurement, core::Vector3f, core::Isometry3f>>& data,
                 bool noise = false);

    inline const SE3PlanarTrajectory& trajectory() const {
      return *trajectory_;
    }

    inline const float freq() const {
      return freq_;
    }

    inline const core::Vector3f& bias_acc() const {
      return ba_;
    }
    inline const core::Vector3f& bias_gyro() const {
      return bg_;
    }
    inline void bias_acc(const core::Vector3f& ba) {
      ba_ = ba;
    }
    inline void bias_gyro(const core::Vector3f& bg) {
      bg_ = bg;
    }

    inline float& std_acc() {
      return std_acc_;
    }
    inline float& std_gyro() {
      return std_gyro_;
    }
    inline float& std_bias_acc() {
      return std_bias_acc_;
    }
    inline float& std_bias_gyro() {
      return std_bias_gyro_;
    }

  protected:
    float freq_; // the frequency at which a new measurement becomes available

    float std_acc_       = 0.00175f;
    float std_gyro_      = 0.00175f;
    float std_bias_acc_  = 0.00167f;
    float std_bias_gyro_ = 0.00167f;

    std::shared_ptr<SE3PlanarTrajectory> trajectory_;
    std::mt19937 rnd_gen_;

    core::Vector3f ba_ = core::Vector3f::Zero();
    core::Vector3f bg_ = core::Vector3f::Zero();
  };

} // namespace test_imu