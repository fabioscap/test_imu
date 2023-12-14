// from KITTI benchmark odometry evaluation

#include <iostream>
#include <limits>
#include <math.h>
#include <stdio.h>
#include <vector>

using namespace std;

#include "common/common.h"
#include <Eigen/Dense>
#include <fstream>

// static parameter
// float lengths[] = {5,10,50,100,150,200,250,300,350,400};
struct errors {
  size_t first_frame;
  float r_err;
  float t_err;
  float len;
  float speed;
  errors(size_t first_frame, float r_err, float t_err, float len, float speed) :
    first_frame(first_frame), r_err(r_err), t_err(t_err), len(len), speed(speed) {
  }
};

vector<std::pair<double, Eigen::Isometry3f>> loadPoses(string file_name) {
  vector<std::pair<double, Eigen::Isometry3f>> poses;
  FILE* fp = fopen(file_name.c_str(), "r");
  if (!fp)
    return poses;
  while (!feof(fp)) {
    Eigen::Isometry3f P;
    P.setIdentity();

    float a, b, c;
    double timestamp;
    if (!fscanf(fp,
                "%lf,%f,%f,%f,%f,%f,%f",
                &timestamp,
                &P.translation().x(),
                &P.translation().y(),
                &P.translation().z(),
                &a,
                &b,
                &c))
      continue;
    P.linear() = test_imu::Rx(a) * test_imu::Ry(b) * test_imu::Rz(c);
    poses.push_back(std::make_pair(timestamp, P));
  }
  fclose(fp);
  return poses;
}

void syncPoses(const vector<std::pair<double, Eigen::Isometry3f>>& p1_stamped,
               const vector<std::pair<double, Eigen::Isometry3f>>& p2_stamped,
               vector<Eigen::Isometry3f>& p1,
               vector<Eigen::Isometry3f>& p2) {
  // downsample the largest file
  const auto& large_stamped = (p1_stamped.size() >= p2_stamped.size()) ? p1_stamped : p2_stamped;
  const auto& small_stamped = (p1_stamped.size() >= p2_stamped.size()) ? p2_stamped : p1_stamped;
  auto& large               = (p1.size() >= p2.size()) ? p1 : p2;
  auto& small               = (p1.size() >= p2.size()) ? p2 : p1;

  size_t j = 0;
  for (size_t i = 0; i < small_stamped.size(); ++i) {
    small.push_back(small_stamped.at(i).second);

    double t = small_stamped.at(i).first;

    while (t > large_stamped.at(j).first) {
      j++;
    }
    // TODO I should check that the timestamps are similar otherwise interpolate the poses
    large.push_back(large_stamped.at(j).second);
  }

  std::cout << p1.size() << "\n" << p2.size() << "\n";
}

vector<float> trajectoryDistances(const vector<Eigen::Isometry3f>& poses) {
  vector<float> dist;
  dist.push_back(0.0f);
  for (size_t i = 1; i < poses.size(); i++) {
    const Eigen::Isometry3f& P1 = poses[i - 1];
    const Eigen::Isometry3f& P2 = poses[i];
    dist.push_back(dist[i - 1] + (P1.translation() - P2.translation()).norm());
  }
  return dist;
}

size_t lastFrameFromSegmentLength(vector<float>& dist, size_t first_frame, float len) {
  for (size_t i = first_frame; i < dist.size(); i++)
    if (dist[i] > dist[first_frame] + len)
      return i;
  return -1;
}

inline float rotationError(const Eigen::Isometry3f& pose_error) {
  const Eigen::Matrix3f R = pose_error.rotation();
  float a                 = R(0, 0);
  float b                 = R(1, 1);
  float c                 = R(2, 2);
  float d                 = 0.5 * (a + b + c - 1.0);
  return acos(max(min(d, 1.0f), -1.0f));
}

inline float translationError(const Eigen::Isometry3f& pose_error) {
  return pose_error.translation().norm();
}

vector<errors> calcSequenceErrors(vector<Eigen::Isometry3f>& poses_gt,
                                  vector<Eigen::Isometry3f>& poses_result,
                                  int len) {
  // error vector
  vector<errors> err;

  // parameters
  size_t step_size = 1;

  // pre-compute distances (from ground truth as reference)
  vector<float> dist = trajectoryDistances(poses_gt);

  // for all start positions do
  for (size_t first_frame = 0; first_frame < poses_gt.size(); first_frame += step_size) {
    // for all segment lengths do
    for (size_t i = 0; i < 1; i++) {
      // compute last frame
      size_t last_frame = lastFrameFromSegmentLength(dist, first_frame, len);

      // continue, if sequence not long enough
      if (last_frame == -1)
        continue;
      // compute rotational and translational errors
      Eigen::Isometry3f pose_delta_gt = (poses_gt[first_frame]).inverse() * poses_gt[last_frame];
      Eigen::Isometry3f pose_delta_result =
        (poses_result[first_frame]).inverse() * poses_result[last_frame];

      Eigen::Isometry3f pose_error = (pose_delta_result).inverse() * pose_delta_gt;
      float r_err                  = rotationError(pose_error);
      float t_err                  = translationError(pose_error);

      // compute speed
      float num_frames = (float) (last_frame - first_frame + 1);
      float speed      = len / (0.1 * num_frames);

      // write to file
      err.push_back(errors(first_frame, r_err / len, t_err / len, len, speed));
    }
  }

  // return error vector
  return err;
}

int main(int argc, char* argv[]) {
  // ground truth and result directories
  if (argc != 4) {
    throw std::runtime_error("compute_rpe len pred.csv gt.csv");
  }

  int len = std::stoi(argv[1]);
  string gt_dir(argv[3]);
  string result_dir(argv[2]);

  // read ground truth and result poses
  std::cout << "loading gt\n";
  vector<std::pair<double, Eigen::Isometry3f>> poses_gt_stamped = loadPoses(gt_dir);
  std::cout << "gt loaded\n";
  std::cout << "loading prediction\n";
  vector<std::pair<double, Eigen::Isometry3f>> poses_result_stamped = loadPoses(result_dir);
  std::cout << "prediction loaded\n";

  std::cout << "sync poses\n";
  vector<Eigen::Isometry3f> poses_gt, poses_result;
  syncPoses(poses_gt_stamped, poses_result_stamped, poses_gt, poses_result);
  std::cout << "done\n";

  // compute sequence errors
  std::cout << "computing errors\n";
  vector<errors> seq_err = calcSequenceErrors(poses_gt, poses_result, len);

  std::string t_err_filename = "t_err_" + std::string(argv[1]) + ".txt";
  std::string r_err_filename = "r_err_" + std::string(argv[1]) + ".txt";
  ofstream t_err(t_err_filename);
  ofstream r_err(r_err_filename);

  for (size_t i = 0; i < seq_err.size(); ++i) {
    const errors& e = seq_err.at(i);
    t_err << e.t_err << "\n";
    r_err << e.r_err << "\n";
  }

  // success
  return 0;
}
