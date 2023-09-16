#include <srrg_geometry/geometry_defs.h>

#include "common/common.h"

#include <ctime>
#include <memory>
#include <random>

// generate an SE(3) trajectory
// and produce imu measurements

namespace test_imu {

  class SE3PlanarTrajectory {
  public:
    SE3PlanarTrajectory(float T) : T_(T), scaling_(2 * M_PI / T_){};

    virtual void sampleTrajectory(float t,
                                  srrg2_core::Vector3f& pos,
                                  srrg2_core::Vector3f& vel,
                                  srrg2_core::Vector3f& acc) const {};

    // The frame is always tangent to the curve with its x-axis,
    // and the z-axis always points up (you can add rotation if needed).
    void getPoseMeasurement(float t, srrg2_core::Isometry3f& pose, IMUMeasurement& measurement);

    inline void setT(float T) {
      T_       = T;
      scaling_ = 2 * M_PI / T_;
    }

    inline float T() const {
      return T_;
    }

  protected:
    float T_;       // The duration of the trajectory
    float scaling_; // linear timing law
  };

  class SE3EightTrajectory : public SE3PlanarTrajectory {
  public:
    SE3EightTrajectory(float T) : SE3PlanarTrajectory(T) {
    }

    void sampleTrajectory(float t,
                          srrg2_core::Vector3f& pos,
                          srrg2_core::Vector3f& vel,
                          srrg2_core::Vector3f& acc) const override;
  };
  class SE3CircleTrajectory : public SE3PlanarTrajectory {
  public:
    SE3CircleTrajectory(float T) : SE3PlanarTrajectory(T) {
    }

    void sampleTrajectory(float t,
                          srrg2_core::Vector3f& pos,
                          srrg2_core::Vector3f& vel,
                          srrg2_core::Vector3f& acc) const override;
  };

  class FakeIMU : public IMUSensor {
  public:
    FakeIMU(std::shared_ptr<SE3PlanarTrajectory> traj_, float freq, int seed = std::time(0)) :
      IMUSensor(freq), trajectory_(traj_), rnd_gen_(seed) {
      srrg2_core::Vector3f pos, vel, acc;
    }

    void generateData(std::vector<std::pair<IMUMeasurement, srrg2_core::Isometry3f>>& data);

    inline const SE3PlanarTrajectory& trajectory() const {
      return *trajectory_;
    }

  protected:
    const float noise_acc_       = 0.00175f;
    const float noise_gyro_      = 0.00175f;
    const float noise_bias_acc_  = 0.00167f;
    const float noise_bias_gyro_ = 0.00167f;

    std::shared_ptr<SE3PlanarTrajectory> trajectory_;
    std::mt19937 rnd_gen_;
  };

} // namespace test_imu