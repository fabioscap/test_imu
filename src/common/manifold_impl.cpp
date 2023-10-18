#include "manifold.h"

namespace test_imu {

  template <typename Manifold>
  ManifoldComp_<Manifold>::ManifoldComp_(const Manifold& m) : manifold_(m) {
  }

  ManifoldSO3::ManifoldSO3() : BaseType(DataType::Identity()) {
  }
  ManifoldSO3::ManifoldSO3(const DataType& data) : BaseType(data) {
  }

  ManifoldSO3 ManifoldSO3::boxplus(const ManifoldSO3::TangentType& dsp) const {
    DataType R = data() * core::geometry3d::expMapSO3(dsp);
    core::fixRotation(R);
    return ManifoldSO3(R);
  }
  ManifoldSO3::TangentType ManifoldSO3::boxminus(const ManifoldSO3& to) const {
    DataType R = data_.transpose() * to.data();
    return core::geometry3d::logMapSO3(R.eval());
  }

  template <int dim>
  Euclidean_<dim> Euclidean_<dim>::boxplus(const Euclidean_<dim>::TangentType& dsp) const {
    return Euclidean_(this->data_ + dsp);
  }
  template <int dim>
  typename Euclidean_<dim>::TangentType Euclidean_<dim>::boxminus(const Euclidean_<dim>& to) const {
    return to.data() - this->data_;
  }
} // namespace test_imu