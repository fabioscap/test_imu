#include <Eigen/Dense>
#include <iostream>

int main() {
  Eigen::Matrix<float, 10, 10> big_matrix = Eigen::Matrix<float, 10, 10>::Zero();

  Eigen::Block<Eigen::Matrix<float, 10, 10>, 3, 3> merda = big_matrix.block<3, 3>(1, 1);

  merda(1, 1) = 10;

  std::cout << merda * Eigen::Vector3f::Ones() << "\n";
  std::cout << big_matrix << "\n" << merda << std::endl;
}