#include "imu_preintegrator/imu_preintegrator.h"

using namespace test_imu;

int main() {
  ImuPreintegrator p;

  ImuMeasurement m;
  float dt = 0.1;

  p.reset();
  p.preintegrate(m, dt);

  return 0;
}
