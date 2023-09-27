#include "synthetic.h"

#include <srrg_geometry/ad.h>
#include <srrg_geometry/geometry3d.h>
#include <srrg_geometry/geometry_defs.h>

#include <cmath>

#include <fstream>
#include <memory>

#include "common/common.h"

using namespace test_imu;

using DVector3f = core::Vector3_<core::DualValuef>;

const float a = 100.0; // amplitude parameter
void SE3EightTrajectory::sampleTrajectory(float t,
                                          core::Vector3f& pos,
                                          core::Vector3f& vel,
                                          core::Vector3f& acc) const {
  // lemniscate cartesian equations
  float scaling = (2 * M_PI / T_);
  float theta   = t * scaling;

  float ctheta  = std::cos(theta);
  float stheta  = std::sin(theta);
  float stheta2 = stheta * stheta;
  float ctheta2 = ctheta * ctheta;

  using namespace core;

  pos << a * Vector3f(ctheta, ctheta * stheta, 0.0) / (stheta2 + (1));
  vel << scaling * a * Vector3f((stheta * (stheta2 - 3)), -(3 * stheta2 - 1), 0) /
           ((stheta2 + 1) * (stheta2 + 1));
  acc << scaling * scaling * a *
           Vector3f((ctheta * (10 * ctheta2 + ctheta2 * ctheta2 - 8)) /
                      ((ctheta2 - 2) * (ctheta2 - 2) * (ctheta2 - 2)),
                    -2 * ctheta * stheta * (3 * ctheta2 + 2) /
                      ((stheta2 + 1) * (stheta2 + 1) * (stheta2 + 1)),
                    0);
}

const float r = 10.0;
void test_imu::SE3CircleTrajectory::sampleTrajectory(float t,
                                                     core::Vector3f& pos,
                                                     core::Vector3f& vel,
                                                     core::Vector3f& acc) const {
  using namespace core;
  float scaling = (2 * M_PI / T_);
  float theta   = t * scaling;

  float ctheta = std::cos(theta);
  float stheta = std::sin(theta);

  pos << r * Vector3f(ctheta, stheta, 0.0);
  vel << scaling * r * Vector3f(-stheta, ctheta, 0.0);
  acc << scaling * scaling * r * Vector3f(-ctheta, -stheta, 0.0);
}

const float l     = 10;
const float theta = 0;
const float v0    = 0.1;
void test_imu::SE3StraightTrajectory::sampleTrajectory(float t,
                                                       srrg2_core::Vector3f& pos,
                                                       srrg2_core::Vector3f& vel,
                                                       srrg2_core::Vector3f& acc) const {
  using namespace core;

  core::Vector3f end(l * std::cos(theta), l * sin(theta), 0);
  float scaling = 1 / T_;
  float s       = t / T_;
  float s2      = s * s;
  float s3      = s2 * s;
  float s4      = s3 * s;
  float s5      = s4 * s;

  pos << (6 * (1 - v0) * s5 - 15 * (1 - v0) * s4 + 10 * (1 - v0) * s3 + v0 * s) * (end);
  vel << (30 * (1 - v0) * s4 - 60 * (1 - v0) * s3 + 30 * (1 - v0) * s2 + v0) * (end) *scaling;
  acc << (120 * (1 - v0) * s3 - 180 * (1 - v0) * s2 + 60 * (1 - v0) * s) * (end) *scaling * scaling;

  // DEBUG version: imu measurements all zero
  // pos << s * (end);
  // vel << scaling * (end);
  // acc << scaling * scaling * 0 * (end);
}

void SE3PlanarTrajectory::getPoseMeasurement(float t,
                                             core::Isometry3f& pose,
                                             ImuMeasurement& measurement) {
  using namespace core;
  pose.setIdentity();

  Vector3f pos, vel, acc;

  sampleTrajectory(t, pos, vel, acc);

  // Pose

  pose.translation() << pos;
  // compute the tangent vector to the curve
  float v_norm = vel.norm();

  if ((std::fabs(v_norm) < 1e-8))
    throw std::runtime_error("velocity appears to be zero ");
  Vector3f tangent = vel / v_norm;

  Vector3f z_vers(0, 0, 1);
  z_vers.normalize();

  Vector3f y_vers = z_vers.cross(tangent);
  y_vers.normalize();

  Eigen::Matrix3f rotation;
  rotation.col(0) << tangent;
  rotation.col(1) << y_vers;
  rotation.col(2) << z_vers;

  pose.linear() = rotation;

  // Measurement
  // linear acceleration in moving frame
  const Matrix3f& R        = pose.rotation();
  measurement.acceleration = R.transpose() * acc;

  // angular velocity in moving frame
  // Rt *Rdot = skew(omega)

  /*   Matrix3f Rdot = Matrix3f::Zero();
    float v_norm2 = v_norm * v_norm;
    Rdot.col(0) << ((Matrix3f::Identity() - vel * vel.transpose() / v_norm2) / v_norm) * acc;
    Rdot.col(1) << -Rdot(1, 0), Rdot(0, 0), 0;
    Matrix3f omega_skew = R.transpose() * Rdot; */

  // expecting only wz != 0
  measurement.angular_vel << 0, 0, (vel(0) * acc(1) - vel(1) * acc(0)) / (v_norm * v_norm);
  measurement.timestamp = t;
}

void test_imu::FakeImu::generateData(std::vector<std::pair<ImuMeasurement, core::Isometry3f>>& data,
                                     bool noise) {
  using namespace core;

  int num  = trajectory_->T() * freq_;
  float dt = 1 / freq_;

  data.resize(num);

  // initialize biases
  Vector3f ba = 1e-8 * Vector3f::Ones();
  Vector3f bg = 1e-8 * Vector3f::Ones();

  std::normal_distribution<double> std_dist(0.0, 1.0);
  std::uniform_real_distribution<double> uni_dist(-1.0, 1.0);

  float t = 0;
  for (int i = 0; i < num; ++i, t += dt) {
    Isometry3f& pose     = data.at(i).second;
    ImuMeasurement& meas = data.at(i).first;

    trajectory_->getPoseMeasurement(t, pose, meas);
    // noise
    // bias is assumed constant between keyframes in 2015 paper (TODO)
    if (noise) {
      ba += noise_bias_acc_ * Vector3f(std_dist(rnd_gen_), std_dist(rnd_gen_), std_dist(rnd_gen_));
      bg += noise_bias_gyro_ * Vector3f(std_dist(rnd_gen_), std_dist(rnd_gen_), std_dist(rnd_gen_));
      meas.acceleration +=
        ba + noise_acc_ * Vector3f(uni_dist(rnd_gen_), uni_dist(rnd_gen_), uni_dist(rnd_gen_));
      meas.angular_vel +=
        bg + noise_gyro_ * Vector3f(uni_dist(rnd_gen_), uni_dist(rnd_gen_), uni_dist(rnd_gen_));
    }
  }
}