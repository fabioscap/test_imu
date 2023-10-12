#include "common.h"
#include <srrg_geometry/geometry3d.h>
namespace test_imu {

  template <int dim>
  using TangentType_ = core::Vector_<float, dim>;

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
    TangentType boxminus(const ThisType&) const {
      throw std::runtime_error("boxminus not defined");
    }

  protected:
    DataType data_;
  };

  class ManifoldSO3 : public ManifoldBase_<core::Matrix3f, 3> {
  public:
    using BaseType = ManifoldBase_<core::Matrix3f, 3>;
    ManifoldSO3() : BaseType(core::Matrix3f::Identity()) {
    }
    ManifoldSO3(const DataType& data) : BaseType(data) {
    }

    ManifoldSO3 boxplus(const TangentType& dsp) const {
      DataType R = data() * core::geometry3d::expMapSO3(dsp);
      core::fixRotation(R);
      return ManifoldSO3(R);
    }
    TangentType boxminus(const ManifoldSO3& to) const {
      DataType R = data_.transpose() * to.data();
      return core::geometry3d::logMapSO3(R.eval());
    }
  };

  template <typename... Manifolds>
  class ManifoldComp_;

  template <typename Manifold>
  class ManifoldComp_<Manifold> {
  public:
    static constexpr int dim = Manifold::dim;

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
    template <size_t N>
    auto& get() {
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