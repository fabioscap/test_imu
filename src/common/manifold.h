#pragma once

#include "common.h"
#include <srrg_geometry/geometry3d.h>

// TODO move as much definitions as possible into impl

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
    ThisType boxplus(const TangentType&) const;

    // convention for boxminus: compute chart around *this
    TangentType boxminus(const ThisType&) const;

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

    ThisType boxplus(const TangentType& dsp) const;

    // convention for boxminus: compute chart around *this
    TangentType boxminus(const ThisType& to) const;

    template <size_t N>
    const Manifold& get() const;

    template <size_t N>
    Manifold& get();

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
    ThisType boxplus(const TangentType& dsp) const;

    // convention for boxminus: compute chart around *this
    TangentType boxminus(const ThisType& to) const;

    template <size_t N>
    auto& get();

    template <size_t N>
    const auto& get() const;

  protected:
    Manifold manifold_;
    ManifoldComp_<Rest...> rest_;
  };

} // namespace test_imu