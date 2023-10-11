#include "common.h"
#include <srrg_geometry/geometry3d.h>
namespace test_imu {

  template <int dim>
  using TangentType_ = core::Vector_<float, dim>;

  template <typename DataType_, int dim_>
  class ManifoldBase {
  public:
    static constexpr int dim = dim_;
    using ThisType           = ManifoldBase<DataType_, dim_>;
    using DataType           = DataType_;
    using TangentType        = TangentType_<dim>;

    inline void setData(const DataType& data) {
      data_ = data;
    }
    inline const DataType& data() const {
      return data_;
    }

    // redefine these two methods
    DataType boxplus(const TangentType&) {
      throw std::runtime_error("boxplus not defined");
    }
    TangentType boxminus(const ThisType&) {
      throw std::runtime_error("boxminus not defined");
    }

  protected:
    DataType data_;
  };

  class ManifoldSO3 : public ManifoldBase<core::Matrix3f, 3> {
  public:
    DataType boxplusSO3(const TangentType& dsp) {
      DataType R = data() * core::geometry3d::expMapSO3(dsp);
      core::fixRotation(R);
      return R;
    }
    TangentType boxminusSO3(const ThisType& to) {
      DataType R = data_.transpose() * to.data();
      return core::geometry3d::logMapSO3(R.eval());
    }
  };

  template <typename... Manifolds>
  class ManifoldComp;

  template <typename Manifold>
  class ManifoldComp<Manifold> {
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
  class ManifoldComp<Manifold, Rest...> {
  public:
    static constexpr int dim = Manifold::dim + ManifoldComp<Rest...>::dim;
    template <size_t N>
    auto& get() {
      if constexpr (N == 0)
        return manifold_;
      else
        return rest_.template get<N - 1>();
    }

  protected:
    Manifold manifold_;
    ManifoldComp<Rest...> rest_;
  };

} // namespace test_imu