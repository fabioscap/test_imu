#include "manifold.h"

namespace test_imu {

  template <typename DataType_, int dim_>
  ManifoldBase_<DataType_, dim_> ManifoldBase_<DataType_, dim_>::boxplus(const TangentType&) const {
    throw std::runtime_error("boxplus not defined");
  }

  template <typename DataType_, int dim_>
  typename ManifoldBase_<DataType_, dim_>::TangentType
  ManifoldBase_<DataType_, dim_>::boxminus(const ThisType&) const {
    throw std::runtime_error("boxminus not defined");
  }

  template <typename Manifold>
  ManifoldComp_<Manifold>::ManifoldComp_(const Manifold& m) : manifold_(m) {
  }

  ManifoldSO3::ManifoldSO3() : BaseType(DataType::Identity()) {
  }
  ManifoldSO3::ManifoldSO3(const DataType& data) : BaseType(data) {
  }

  ManifoldSO3 ManifoldSO3::boxplus(const TangentType& dsp) const {
    DataType R = data() * core::geometry3d::expMapSO3(dsp);
    core::fixRotation(R);
    return ManifoldSO3(R);
  }
  ManifoldSO3::TangentType ManifoldSO3::boxminus(const ManifoldSO3& to) const {
    DataType R = data_.transpose() * to.data();
    return core::geometry3d::logMapSO3(R.eval());
  }

  ManifoldSO2::ManifoldSO2() : BaseType(0.0) {
  }
  ManifoldSO2::ManifoldSO2(const DataType& data) : BaseType(data) {
  }

  ManifoldSO2 ManifoldSO2::boxplus(const TangentType& dsp) const {
    Scalar s = std::sin(data() + dsp(0));
    Scalar c = std::cos(data() + dsp(0));
    return ManifoldSO2(std::atan2(s, c));
  }

  ManifoldSO2::TangentType ManifoldSO2::boxminus(const ManifoldSO2& to) const {
    Scalar s = std::sin(to.data() - data());
    Scalar c = std::cos(to.data() - data());
    return TangentType(std::atan2(s, c));
  }

  template <int dim>
  Euclidean_<dim> Euclidean_<dim>::boxplus(const TangentType& dsp) const {
    return Euclidean_(this->data_ + dsp);
  }
  template <int dim>
  typename Euclidean_<dim>::TangentType Euclidean_<dim>::boxminus(const Euclidean_<dim>& to) const {
    return to.data() - this->data_;
  }

  template <typename Manifold>
  ManifoldComp_<Manifold> ManifoldComp_<Manifold>::boxplus(const TangentType& dsp) const {
    return ManifoldComp_<Manifold>(manifold_.boxplus(dsp));
  }

  template <typename Manifold>
  template <size_t N>
  const Manifold& ManifoldComp_<Manifold>::get() const {
    static_assert(N == 0, "Index out of range");
    return manifold_;
  }

  template <typename Manifold>
  template <size_t N>
  Manifold& ManifoldComp_<Manifold>::get() {
    static_assert(N == 0, "Index out of range");
    return manifold_;
  }

  template <typename Manifold>
  typename ManifoldComp_<Manifold>::TangentType
  ManifoldComp_<Manifold>::boxminus(const ThisType& to) const {
    return manifold_.boxminus(to.get<0>());
  }

  template <typename Manifold, typename... Rest>
  ManifoldComp_<Manifold, Rest...>
  ManifoldComp_<Manifold, Rest...>::boxplus(const TangentType& dsp) const {
    ManifoldComp_<Manifold, Rest...> out;
    // the first
    out.get<0>() = manifold_.boxplus(dsp.head(manifold_.dim));
    // the rest
    out.rest_ = rest_.boxplus(dsp.tail(dim - manifold_.dim));
    return out;
  }

  template <typename Manifold, typename... Rest>
  typename ManifoldComp_<Manifold, Rest...>::TangentType
  ManifoldComp_<Manifold, Rest...>::boxminus(const ThisType& to) const {
    ManifoldComp_<Manifold, Rest...>::TangentType out;
    // the first
    out.head(manifold_.dim) = manifold_.boxminus(to.get<0>());
    // the rest
    out.tail(dim - manifold_.dim) = rest_.boxminus(to.rest_);
    return out;
  }

  template <typename Manifold, typename... Rest>
  template <size_t N>
  auto& ManifoldComp_<Manifold, Rest...>::get() {
    if constexpr (N == 0)
      return manifold_;
    else
      return rest_.template get<N - 1>();
  }

  template <typename Manifold, typename... Rest>
  template <size_t N>
  const auto& ManifoldComp_<Manifold, Rest...>::get() const {
    if constexpr (N == 0)
      return manifold_;
    else
      return rest_.template get<N - 1>();
  }
} // namespace test_imu