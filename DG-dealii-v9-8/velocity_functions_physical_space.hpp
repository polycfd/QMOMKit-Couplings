#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

/**
 * @brief Zero velocity field. Just a wrapper around ConstantFunction
 */
template <int dim>
using NoFlow = ZeroTensorFunction<1, dim>;

/**
 * @brief Constant flow along X direction
 */
template <int dim>
class ConstantXFlow : public TensorFunction<1, dim>
{
public:
  /**
   * @brief Constructor.
   *
   */
  ConstantXFlow()
    : TensorFunction<1, dim>()
  {}

  /**
   * @brief Return the velocity
   *
   * @param p         Spatial point
   * @param component Component of the velocity field.
   */
  virtual Tensor<1, dim>
  value(const Point<dim> & /*p*/) const override
  {
    Tensor<1, dim> vel;
    vel[0] = 1.;
    return vel;
  }
};

/**
 * Single-cell Taylor–Green vortex on [-1,1]^2, divergence-free,
 * with zero normal velocity on the boundary (free-slip).
 *
 *   u(x,y) = -A * cos(pi*x/2) * sin(pi*y/2)
 *   v(x,y) =  A * sin(pi*x/2) * cos(pi*y/2)
 */
template <int dim>
class TaylorGreenVortex : public TensorFunction<1, dim>
{
public:
  /**
   * @brief Constructor.
   *
   */
  TaylorGreenVortex()
    : TensorFunction<1, dim>()
  {}

  /**
   * @brief Return the velocity
   *
   * @param p         Spatial point
   */
  virtual Tensor<1, dim>
  value(const Point<dim> &p) const override
  {
    Tensor<1, dim> vel;

    // The vortex only exists in 2D. The other instantiations still have to
    // compile, since the problem is picked at run time, but they can never
    // be reached.
    if constexpr (dim == 2)
      {
        const double k  = 0.5 * numbers::PI;
        const double cx = std::cos(k * p[0]);
        const double sx = std::sin(k * p[0]);
        const double cy = std::cos(k * p[1]);
        const double sy = std::sin(k * p[1]);

        vel[0] = -omega * cx * sy;
        vel[1] = omega * sx * cy;
      }
    else
      {
        AssertThrow(false,
                    ExcMessage("The Taylor-Green vortex is only defined in "
                               "two dimensions."));
      }
    return vel;
  }

  const double omega = 1;
};