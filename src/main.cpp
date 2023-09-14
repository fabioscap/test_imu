#include <Eigen/Dense>
#include <srrg_geometry/ad.h>

#include <srrg_solver/solver_core/factor_graph.h>

#include <iostream>

int main() {
  Eigen::Matrix3f b = Eigen::Matrix3f::Zero();
  srrg2_core::DualValuef a(5);

  srrg2_solver::FactorGraph* penis = new srrg2_solver::FactorGraph();

  std::cout << "ciaoo" << std::endl;

  return 0;
}