#pragma once

#include "common.h"
#include <srrg_geometry/geometry3d.h>
namespace test_imu {

  template <int dim>
  using TangentType_ = core::Vector_<Scalar, dim>;

  template <typename DataType_, int dim_>
  class ManifoldBase_ {
  public:
    static constexpr int dim = dim_;
    using ThisType           = ManifoldBase_<DataType_, dim_>;
    using DataType           = DataType_;
    using TangentType        = TangentType_<dim>;
    ManifoldBase_()          = default;
    ManifoldBase_(const DataType& data) : data_(data) {
    }

    inline void setData(const DataType& data) {
      data_ = data;
    }
    inline const DataType& data() const {
      return data_;
    }

    // redefine these two methods
    ThisType boxplus(const TangentType&) const {
      throw std::runtime_error("boxplus not defined");
    }

    // convention for boxminus: compute chart around *this
    TangentType boxminus(const ThisType&) const {
      throw std::runtime_error("boxminus not defined");
    }

  protected:
    DataType data_;
  };

  class ManifoldSO3 : public ManifoldBase_<core::Matrix3_<Scalar>, 3> {
  public:
    using BaseType = ManifoldBase_<core::Matrix3_<Scalar>, 3>;
    ManifoldSO3();
    ManifoldSO3(const DataType& data);
    ManifoldSO3 boxplus(const TangentType& dsp) const;
    TangentType boxminus(const ManifoldSO3& to) const;
  };

  // specialization for Euclidean domains
  template <int dim>
  class Euclidean_ : public ManifoldBase_<TangentType_<dim>, dim> {
  public:
    using BaseType    = ManifoldBase_<TangentType_<dim>, dim>;
    using DataType    = typename BaseType::DataType;
    using TangentType = typename BaseType::TangentType;
    Euclidean_() : BaseType(DataType::Zero()) {
    }
    Euclidean_(const DataType& data) : BaseType(data) {
    }
    Euclidean_ boxplus(const TangentType& dsp) const;

    TangentType boxminus(const Euclidean_& to) const;
  };

  template <typename... Manifolds>
  class ManifoldComp_;

  template <typename Manifold>
  class ManifoldComp_<Manifold> {
  public:
    static constexpr int dim = Manifold::dim;
    using ThisType           = ManifoldComp_<Manifold>;
    using TangentType        = typename Manifold::TangentType;

    ManifoldComp_() = default;

    ManifoldComp_(const Manifold& m);

    ThisType boxplus(const TangentType& dsp) const {
      return ThisType(manifold_.boxplus(dsp));
    }

    // convention for boxminus: compute chart around *this
    TangentType boxminus(const ThisType& to) const {
      return manifold_.boxminus(to.get<0>());
    }

    template <size_t N>
    const Manifold& get() const {
      static_assert(N == 0, "Index out of range");
      return manifold_;
    }

    template <size_t N>
    Manifold& get() {
      static_assert(N == 0, "Index out of range");
      return manifold_;
    }

  protected:
    Manifold manifold_;
  };

  template <typename Manifold, typename... Rest>
  class ManifoldComp_<Manifold, Rest...> {
  public:
    static constexpr int dim = Manifold::dim + ManifoldComp_<Rest...>::dim;
    using ThisType           = ManifoldComp_<Manifold, Rest...>;
    using TangentType        = TangentType_<dim>;

    // I hope compiler does return value optimization
    ThisType boxplus(const TangentType& dsp) const {
      ThisType out;
      // the first
      out.get<0>() = manifold_.boxplus(dsp.head(manifold_.dim));
      // the rest
      out.rest_ = rest_.boxplus(dsp.tail(dim - manifold_.dim));
      return out;
    }

    // convention for boxminus: compute chart around *this
    TangentType boxminus(const ThisType& to) const {
      TangentType out;
      // the first
      out.head(manifold_.dim) = manifold_.boxminus(to.get<0>());
      // the rest
      out.tail(dim - manifold_.dim) = rest_.boxminus(to.rest_);
      return out;
    }

    template <size_t N>
    auto& get() {
      if constexpr (N == 0)
        return manifold_;
      else
        return rest_.template get<N - 1>();
    }

    template <size_t N>
    const auto& get() const {
      if constexpr (N == 0)
        return manifold_;
      else
        return rest_.template get<N - 1>();
    }

  protected:
    Manifold manifold_;
    ManifoldComp_<Rest...> rest_;
  };

} // namespace test_imu