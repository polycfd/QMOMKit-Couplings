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

using namespace dealii;

/**
 * @brief Time-dependent analytical raw moments of a Gaussian distribution
 * that is constant in physical space.
 *
 * This class inherits from deal.II's Function<dim> and provides the
 * analytical solution for the raw moments of a normal distribution
 * \f$X \sim \mathcal{N}(\mu, \sigma^2)\f$. Each vector component @p k
 * returns the \f$k\f$-th raw moment scaled by an exponential in time:
 * \f[
 *   M_k(t) = \mathrm{E}[X^{k}] \, e^{k\,t}
 * \f]
 * where the raw moment of order \f$n = k\f$ is computed from the
 * general formula for a Gaussian distribution:
 * \f[
 *   \mathrm{E}[X^n] = \sum_{j=0}^{\lfloor n/2 \rfloor}
 *     \binom{n}{2j} \, (2j-1)!! \, \sigma^{2j} \, \mu^{n-2j}
 * \f]
 * This supports an arbitrary number of moments.
 *
 * At \f$t = 0\f$ this reduces to the static Gaussian raw moments,
 * making the class suitable both as an initial condition and as an
 * analytical reference solution. The time is obtained from the base
 * class via Function::get_time().
 *
 * The returned value is spatially uniform (independent of the point
 * @p p), since the moment field is homogeneous over the mesh.
 *
 * @tparam dim Spatial dimension of the mesh.
 */
template <int dim>
class ConstantPhysicalGaussianPhase : public Function<dim>
{
public:
  /**
   * @brief Constructor
   *
   * @param n_moments Number of moment equations
   * @param mu  Mean of the Gaussian distribution
   * @param sigma Standard deviation of the Gaussian distribution
   */
  ConstantPhysicalGaussianPhase(const unsigned int n_moments,
                                const double       mu,
                                const double       sigma)
    : Function<dim>(n_moments)
    , mu(mu)
    , sigma(sigma)
  {}

  /**
   * @brief Return the time-dependent (component+1)-th raw moment.
   *
   * @param p Spatial point (unused — the field is uniform).
   * @param component Index selecting which moment to return.
   * @return \f$\mathrm{E}[X^{k}] \, e^{k\,t}\f$ where
   *         \f$k\f$ = @p component and \f$t\f$ is the current time.
   */
  virtual double
  value(const Point<dim> & /*p*/,
        const unsigned int component = 0) const override
  {
    const double       t = this->get_time();
    const unsigned int n = component;

    // Compute E[X^n] for X ~ N(mu_phase, sigma_phase^2) using:
    // E[X^n] = sum_{k=0}^{floor(n/2)} C(n,2k) * (2k-1)!! * sigma_phase^{2k}
    //                                                       * mu_phase^{n-2k}
    double moment = 0.0;
    for (unsigned int k = 0; k <= n / 2; ++k)
      {
        // Binomial coefficient C(n, 2k)
        double binom = 1.0;
        for (unsigned int j = 0; j < 2 * k; ++j)
          binom *= static_cast<double>(n - j) / static_cast<double>(j + 1);

        // Double factorial (2k-1)!! with (-1)!! = 1 by convention
        double dfact = 1.0;
        for (unsigned int j = 1; j < 2 * k; j += 2)
          dfact *= j;

        moment += binom * dfact * std::pow(sigma, 2.0 * k) *
                  std::pow(mu, static_cast<double>(n - 2 * k));
      }

    return moment * std::exp(static_cast<double>(n) * t);
  }

  /// Mean of the Gaussian distribution (phase space)
  const double mu;
  /// Standard deviation of the Gaussian distribution (phase space)
  const double sigma;
};

/**
 * @brief Time-dependent analytical raw moments of a Gaussian distribution that
 * also varies like a gaussian in physical space.
 *
 * This class inherits from deal.II's Function<dim> and provides the
 * analytical solution for the raw moments of a normal distribution
 * \f$X \sim \mathcal{N}(\mu, \sigma^2)\f$. Each vector component @p k
 * returns the \f$k\f$-th raw moment scaled by an exponential in time:
 * \f[
 *   M_k(x,t) = \mathrm{E}[X^{k}] \, e^{k\,t}
 * \f]
 * where the raw moment of order \f$n = k\f$ is computed from the
 * general formula for a Gaussian distribution:
 * \f[
 *   \mathrm{E}[X^n] = \sum_{j=0}^{\lfloor n/2 \rfloor}
 *     \binom{n}{2j} \, (2j-1)!! \, \sigma^{2j} \, \mu^{n-2j}
 * \f]
 * This supports an arbitrary number of moments.
 *
 * At \f$t = 0\f$ this reduces to the static Gaussian raw moments,
 * making the class suitable both as an initial condition and as an
 * analytical reference solution. The time is obtained from the base
 * class via Function::get_time().
 *
 * @tparam dim Spatial dimension of the mesh.
 */
template <int dim>
class GaussianPhysicalGaussianPhase : public Function<dim>
{
public:
  /**
   * @brief Constructor
   *
   * @param n_moments Number of moment equations
   * @param mu_phase Mean of the Gaussian distribution
   * @param sigma_phase Standard deviation of the Gaussian distribution
   * @param sigma_physical Standard deviation of the Gaussian distribution in
   * physical space
   */
  GaussianPhysicalGaussianPhase(const unsigned int n_moments,
                                const double       mu_phase,
                                const double       sigma_phase,
                                const double       sigma_physical)
    : Function<dim>(n_moments)
    , mu_phase(mu_phase)
    , sigma_phase(sigma_phase)
    , sigma_physical(sigma_physical)
  {}

  /**
   * @brief Return the time-dependent (component+1)-th raw moment.
   *
   * @param p         Spatial point (unused — the field is uniform).
   * @param component Index selecting which moment to return.
   * @return \f$\mathrm{E}[X^{k}] \, e^{a*k\,t}\f$ where
   *         \f$k\f$ = @p component and \f$t\f$ is the current time.
   */
  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    const double       t = this->get_time();
    const unsigned int n = component;

    // Compute E[X^n] for X ~ N(mu_phase, sigma_phase^2) using:
    // E[X^n] = sum_{k=0}^{floor(n/2)} C(n,2k) * (2k-1)!! * sigma_phase^{2k}
    //                                                       * mu_phase^{n-2k}
    double moment = 0.0;
    for (unsigned int k = 0; k <= n / 2; ++k)
      {
        // Binomial coefficient C(n, 2k)
        double binom = 1.0;
        for (unsigned int j = 0; j < 2 * k; ++j)
          binom *= static_cast<double>(n - j) / static_cast<double>(j + 1);

        // Double factorial (2k-1)!! with (-1)!! = 1 by convention
        double dfact = 1.0;
        for (unsigned int j = 1; j < 2 * k; j += 2)
          dfact *= j;

        moment += binom * dfact * std::pow(sigma_phase, 2.0 * k) *
                  std::pow(mu_phase, static_cast<double>(n - 2 * k));
      }

    Point<dim> shifted_p = p;
    shifted_p[0] -= t; // shift by velocity_x * t

    // Wrap into [-1, 1]
    const double L = 2.0;
    shifted_p[0]   = shifted_p[0] - L * std::floor((shifted_p[0] + 1.0) / L);

    return moment * std::exp(static_cast<double>(n) * t) *
           std::exp(-shifted_p.norm_square() /
                    (2.0 * sigma_physical * sigma_physical));
  }

  /// Mean of the Gaussian distribution (phase space)
  const double mu_phase;
  /// Standard deviation of the Gaussian distribution (phase space)
  const double sigma_phase;
  /// Standard deviation of the Gaussian distribution (physical space)
  const double sigma_physical;
};

/**
 * @brief Time-dependent analytical raw moments of a Gaussian distribution that
 * also varies like a gaussian in physical space, and with a constant phase
 * space velocity.
 *
 * This class inherits from deal.II's Function<dim> and provides the
 * analytical solution for the raw moments of a normal distribution
 * \f$X \sim \mathcal{N}(\mu, \sigma^2)\f$. Each vector component @p k
 * returns the \f$k\f$-th raw moment scaled by an exponential in time:
 * \f[
 *   M_k(x,t) = \mathrm{E}[X^{k}]
 * \f]
 * where the raw moment of order \f$n = k\f$ is computed from the
 * general formula for a Gaussian distribution:
 * \f[
 *   \mathrm{E}[X^n] = \sum_{j=0}^{\lfloor n/2 \rfloor}
 *     \binom{n}{2j} \, (2j-1)!! \, \sigma^{2j} \, \mu^{n-2j}
 * \f]
 * This supports an arbitrary number of moments.
 *
 * At \f$t = 0\f$ this reduces to the static Gaussian raw moments,
 * making the class suitable both as an initial condition and as an
 * analytical reference solution. The time is obtained from the base
 * class via Function::get_time().
 *
 * @tparam dim Spatial dimension of the mesh.
 */
template <int dim>
class GaussianPhysicalGaussianPhaseConstant : public Function<dim>
{
public:
  /**
   * @brief Constructor
   *
   * @param n_moments Number of moment equations
   * @param mu_phase Mean of the Gaussian distribution
   * @param sigma_phase Standard deviation of the Gaussian distribution
   * @param sigma_physical Standard deviation of the Gaussian distribution in
   * physical space
   */
  GaussianPhysicalGaussianPhaseConstant(const unsigned int n_moments,
                                        const double       mu_phase,
                                        const double       sigma_phase,
                                        const double       sigma_physical)
    : Function<dim>(n_moments)
    , mu_phase(mu_phase)
    , sigma_phase(sigma_phase)
    , sigma_physical(sigma_physical)
  {}

  /**
   * @brief Return the time-dependent (component+1)-th raw moment.
   *
   * @param p         Spatial point (unused — the field is uniform).
   * @param component Index selecting which moment to return.
   * @return \f$\mathrm{E}[X^{k}] \f$ where
   *         \f$k\f$ = @p component and \f$t\f$ is the current time.
   */
  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    const double       t = this->get_time();
    const unsigned int n = component;

    // Compute E[X^n] for X ~ N(mu_phase, sigma_phase^2) using:
    // E[X^n] = sum_{k=0}^{floor(n/2)} C(n,2k) * (2k-1)!! * sigma_phase^{2k}
    //                                                       * mu_phase^{n-2k}

    double phase_space_growth_velocity = 0.1;
    double mean = mu_phase + phase_space_growth_velocity * t;

    double moment = 0.0;
    for (unsigned int k = 0; k <= n / 2; ++k)
      {
        // Binomial coefficient C(n, 2k)
        double binom = 1.0;
        for (unsigned int j = 0; j < 2 * k; ++j)
          binom *= static_cast<double>(n - j) / static_cast<double>(j + 1);

        // Double factorial (2k-1)!! with (-1)!! = 1 by convention
        double dfact = 1.0;
        for (unsigned int j = 1; j < 2 * k; j += 2)
          dfact *= j;

        moment += binom * dfact * std::pow(sigma_phase, 2.0 * k) *
                  std::pow(mean, static_cast<double>(n - 2 * k));
      }

    Point<dim> shifted_p = p;
    shifted_p[0] -= t; // shift by velocity_x * t

    // Wrap into [-1, 1]
    const double L = 2.0;
    shifted_p[0]   = shifted_p[0] - L * std::floor((shifted_p[0] + 1.0) / L);

    return moment * std::exp(-shifted_p.norm_square() /
                             (2.0 * sigma_physical * sigma_physical));
  }

  /// Mean of the Gaussian distribution (phase space)
  const double mu_phase;
  /// Standard deviation of the Gaussian distribution (phase space)
  const double sigma_phase;
  /// Standard deviation of the Gaussian distribution (physical space)
  const double sigma_physical;
};

template <int dim>
class GaussianPhysicalConstantPhase : public Function<dim>
{
public:
  /**
   * @brief Constructor.
   *
   * @param n_moments Number of moment equations
   * @param mu_phase Mean of the Gaussian distribution
   * @param sigma_phase Standard deviation of the Gaussian distribution
   * @param sigma_physical Standard deviation of the Gaussian distribution
   * in physical space
   */
  GaussianPhysicalConstantPhase(const unsigned int n_moments,
                                const double       mu_phase,
                                const double       sigma_phase,
                                const double       sigma_physical)
    : Function<dim>(n_moments)
    , mu_phase(mu_phase)
    , sigma_phase(sigma_phase)
    , sigma_physical(sigma_physical)
  {}

  /**
   * @brief Return the time-dependent (component+1)-th raw moment.
   *
   * @param p         Spatial point (unused — the field is uniform).
   * @param component Index selecting which moment to return.
   * @return \f$\mathrm{E}[X^{k}] \, e^{k\,t}\f$ where
   *         \f$k\f$ = @p component and \f$t\f$ is the current time.
   */
  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    const double       t = this->get_time();
    const unsigned int n = component;

    // Compute E[X^n] for X ~ N(mu_phase, sigma_phase^2) using:
    // E[X^n] = sum_{k=0}^{floor(n/2)} C(n,2k) * (2k-1)!! * sigma_phase^{2k}
    //                                                       * mu_phase^{n-2k}
    double moment = 0.0;
    for (unsigned int k = 0; k <= n / 2; ++k)
      {
        // Binomial coefficient C(n, 2k)
        double binom = 1.0;
        for (unsigned int j = 0; j < 2 * k; ++j)
          binom *= static_cast<double>(n - j) / static_cast<double>(j + 1);

        // Double factorial (2k-1)!! with (-1)!! = 1 by convention
        double dfact = 1.0;
        for (unsigned int j = 1; j < 2 * k; j += 2)
          dfact *= j;

        moment += binom * dfact * std::pow(sigma_phase, 2.0 * k) *
                  std::pow(mu_phase, static_cast<double>(n - 2 * k));
      }

    Point<dim> shifted_p = p;
    shifted_p[0] -= t; // shift by velocity_x * t

    // Wrap into [-1, 1]
    const double L = 2.0;
    shifted_p[0]   = shifted_p[0] - L * std::floor((shifted_p[0] + 1.0) / L);

    return moment * std::exp(-shifted_p.norm_square() /
                             (2.0 * sigma_physical * sigma_physical));
  }

  /// Mean of the Gaussian distribution (phase space)
  const double mu_phase;
  /// Standard deviation of the Gaussian distribution (phase space)
  const double sigma_phase;
  /// Standard deviation of the Gaussian distribution (physical space)
  const double sigma_physical;
};

/**
 * @brief Exact raw moments of a Gaussian blob stirred by an arbitrary
 * steady, divergence-free velocity field, with linear growth in phase
 * space.
 *
 * The moments transported by this code obey
 * \f[
 *   \frac{\partial M_k}{\partial t}
 *     + \nabla_x \cdot \left( \mathbf{u} M_k \right) = k \, M_k ,
 * \f]
 * the source term being closed exactly by \f$G(\xi) = \xi\f$. For a
 * divergence-free \f$\mathbf{u}\f$ this reduces to the advection form
 * \f$\partial_t M_k + \mathbf{u} \cdot \nabla M_k = k M_k\f$, whose
 * characteristics are simply the streamlines. Along a characteristic the
 * moment obeys \f$\mathrm{d}M_k/\mathrm{d}t = k M_k\f$, so the exact
 * solution is the initial field evaluated at the *foot* of the
 * characteristic, amplified by the phase-space growth:
 * \f[
 *   M_k(\mathbf{x}, t)
 *     = M_k\!\left( \mathbf{X}_0(\mathbf{x}, t),\, 0 \right) e^{k t} ,
 * \f]
 * where \f$\mathbf{X}_0(\mathbf{x},t)\f$ is the position at time 0 of the
 * particle located at \f$\mathbf{x}\f$ at time \f$t\f$.
 *
 * The initial field is a Gaussian blob in physical space carrying the raw
 * moments of \f$X \sim \mathcal{N}(\mu, \sigma^2)\f$:
 * \f[
 *   M_k(\mathbf{x}, 0) = \mathrm{E}[X^{k}] \,
 *     \exp\!\left( -\frac{\|\mathbf{x} - \mathbf{x}_c\|^2}
 *                        {2 \sigma_{space}^2} \right) ,
 * \f]
 * \f[
 *   \mathrm{E}[X^n] = \sum_{j=0}^{\lfloor n/2 \rfloor}
 *     \binom{n}{2j} \, (2j-1)!! \, \sigma^{2j} \, \mu^{n-2j} .
 * \f]
 *
 * The Taylor-Green vortex has no elementary closed form for its
 * trajectories (they reduce to incomplete elliptic integrals), so
 * \f$\mathbf{X}_0\f$ is obtained by integrating
 * \f$\mathrm{d}\mathbf{y}/\mathrm{d}s = -\mathbf{u}(\mathbf{y})\f$ from
 * \f$\mathbf{y}(0) = \mathbf{x}\f$ over a pseudo-time \f$s \in [0, t]\f$
 * with classical RK4. This is legitimate because the velocity field is
 * steady: the backward flow map depends only on the elapsed time. The
 * substep is capped at #max_substep, which keeps the trajectory error
 * several orders of magnitude below the discretisation error of the DG
 * scheme being measured.
 *
 * Since the blob shears as it is stirred, it does *not* remain Gaussian —
 * only the composition with the backward flow map is exact.
 *
 * @tparam dim Spatial dimension of the mesh.
 */
template <int dim>
class SwirlGaussianPhysicalGaussianPhase : public Function<dim>
{
public:
  /**
   * @brief Constructor.
   *
   * @param n_moments   Number of moment equations (= number of vector
   * components of the Function).
   * @param mu_phase    Mean of the Gaussian distribution.
   * @param sigma_phase Standard deviation of the Gaussian distribution.
   * @param sigma_physical Width of the Gaussian blob in physical space.
   * @param center      Initial location of the blob.
   * @param velocity    Steady, divergence-free velocity field the blob is
   * stirred by. It must be the very field used to advect the solution,
   * otherwise the reference is not a solution of the discretised problem.
   */
  SwirlGaussianPhysicalGaussianPhase(
    const unsigned int                                  n_moments,
    const double                                        mu_phase,
    const double                                        sigma_phase,
    const double                                        sigma_physical,
    const Point<dim>                                    center,
    const std::shared_ptr<const TensorFunction<1, dim>> velocity)
    : Function<dim>(n_moments)
    , mu_phase(mu_phase)
    , sigma_phase(sigma_phase)
    , sigma_physical(sigma_physical)
    , center(center)
    , velocity(velocity)
  {}

  /**
   * @brief Return the exact (component)-th raw moment at @p p.
   *
   * @param p         Spatial point.
   * @param component Index selecting which moment to return.
   * @return \f$M_k(\mathbf{X}_0(p,t), 0) \, e^{k t}\f$ where
   *         \f$k\f$ = @p component and \f$t\f$ is the current time.
   */
  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    const double       t = this->get_time();
    const unsigned int n = component;

    // Compute E[X^n] for X ~ N(mu_phase, sigma_phase^2) using:
    // E[X^n] = sum_{k=0}^{floor(n/2)} C(n,2k) * (2k-1)!! * sigma_phase^{2k}
    //                                                       * mu_phase^{n-2k}
    double moment = 0.0;
    for (unsigned int k = 0; k <= n / 2; ++k)
      {
        // Binomial coefficient C(n, 2k)
        double binom = 1.0;
        for (unsigned int j = 0; j < 2 * k; ++j)
          binom *= static_cast<double>(n - j) / static_cast<double>(j + 1);

        // Double factorial (2k-1)!! with (-1)!! = 1 by convention
        double dfact = 1.0;
        for (unsigned int j = 1; j < 2 * k; j += 2)
          dfact *= j;

        moment += binom * dfact * std::pow(sigma_phase, 2.0 * k) *
                  std::pow(mu_phase, static_cast<double>(n - 2 * k));
      }

    // Follow the characteristic back to where this material point started,
    // and read the initial blob there.
    const Tensor<1, dim> shifted_p = trace_back_to_origin(p, t) - center;

    return moment * std::exp(static_cast<double>(n) * t) *
           std::exp(-shifted_p.norm_square() /
                    (2.0 * sigma_physical * sigma_physical));
  }

  /// Mean of the Gaussian distribution (phase space)
  const double mu_phase;
  /// Standard deviation of the Gaussian distribution (phase space)
  const double sigma_phase;
  /// Standard deviation of the Gaussian distribution (physical space)
  const double sigma_physical;
  /// Initial space location of the Gaussian
  const Point<dim> center;

private:
  /**
   * @brief Integrate the characteristic backwards from @p p over a time
   * @p t, i.e. evaluate the backward flow map
   * \f$\mathbf{X}_0(\mathbf{x}, t)\f$.
   *
   * Classical RK4 applied to \f$\mathrm{d}\mathbf{y}/\mathrm{d}s =
   * -\mathbf{u}(\mathbf{y})\f$. The number of substeps grows linearly with
   * @p t so that the substep never exceeds #max_substep.
   *
   * The velocity fields used here have zero normal component on the
   * boundary of the domain, so the trajectory can never leave it and no
   * clipping is required.
   */
  Point<dim>
  trace_back_to_origin(const Point<dim> &p, const double t) const
  {
    if (t <= 0.)
      return p;

    const unsigned int n_substeps =
      static_cast<unsigned int>(std::ceil(t / max_substep));
    const double h = t / n_substeps;

    Point<dim> y = p;
    for (unsigned int s = 0; s < n_substeps; ++s)
      {
        const Tensor<1, dim> k1 = -velocity->value(y);
        const Tensor<1, dim> k2 = -velocity->value(y + (0.5 * h) * k1);
        const Tensor<1, dim> k3 = -velocity->value(y + (0.5 * h) * k2);
        const Tensor<1, dim> k4 = -velocity->value(y + h * k3);

        y += (h / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);
      }
    return y;
  }

  /// Largest RK4 substep used when tracing a characteristic. The flow
  /// turns over on a time scale of order 1, so this resolves it with ~20
  /// steps per revolution and leaves the trajectory error near 1e-6.
  static constexpr double max_substep = 0.05;

  /// Velocity field defining the characteristics.
  const std::shared_ptr<const TensorFunction<1, dim>> velocity;
};

/**
 * @brief Time-dependent analytical raw moments of a Gaussian distribution
 * that is constant in physical space and includes a source term.
 *
 * This class inherits from deal.II's Function<dim> and provides the
 * analytical solution for the raw moments of a normal distribution
 * \f$X \sim \mathcal{N}(\mu, \sigma^2)\f$. Each vector component @p k
 * returns the \f$k\f$-th raw moment scaled by an exponential in time:
 * \f[
 *   M_k(t) = \mathrm{E}[X^{k}] \, e^{(k+alpha)\,t}
 * \f]
 * where the raw moment of order \f$n = k\f$ is computed from the
 * general formula for a Gaussian distribution:
 * \f[
 *   \mathrm{E}[X^n] = \sum_{j=0}^{\lfloor n/2 \rfloor}
 *     \binom{n}{2j} \, (2j-1)!! \, \sigma^{2j} \, \mu^{n-2j}
 * \f]
 * This supports an arbitrary number of moments.
 *
 * \f$ \alpha \f$ is the multiplication rate in the source term
 *
 * At \f$t = 0\f$ this reduces to the static Gaussian raw moments,
 * making the class suitable both as an initial condition and as an
 * analytical reference solution. The time is obtained from the base
 * class via Function::get_time().
 *
 * The returned value is spatially uniform (independent of the point
 * @p p), since the moment field is homogeneous over the mesh.
 *
 * @tparam dim Spatial dimension of the mesh.
 */
template <int dim>
class ConstantPhysicalGaussianPhaseWithSource : public Function<dim>
{
public:
  /**
   * @brief Constructor
   *
   * @param n_moments Number of moment equations
   * @param mu  Mean of the Gaussian distribution
   * @param sigma Standard deviation of the Gaussian distribution
   * @param alpha Multiplication rate
   */
  ConstantPhysicalGaussianPhaseWithSource(const unsigned int n_moments,
                                          const double       mu,
                                          const double       sigma,
                                          const double       alpha)
    : Function<dim>(n_moments)
    , mu(mu)
    , sigma(sigma)
    , alpha(alpha)
  {}

  /**
   * @brief Return the time-dependent (component+1)-th raw moment.
   *
   * @param p Spatial point (unused — the field is uniform).
   * @param component Index selecting which moment to return.
   * @return \f$\mathrm{E}[X^{k}] \, e^{k\,t}\f$ where
   *         \f$k\f$ = @p component and \f$t\f$ is the current time.
   */
  virtual double
  value(const Point<dim> & /*p*/,
        const unsigned int component = 0) const override
  {
    const double       t = this->get_time();
    const unsigned int n = component;

    // Compute E[X^n] for X ~ N(mu_phase, sigma_phase^2) using:
    // E[X^n] = sum_{k=0}^{floor(n/2)} C(n,2k) * (2k-1)!! * sigma_phase^{2k}
    //                                                       * mu_phase^{n-2k}
    double moment = 0.0;
    for (unsigned int k = 0; k <= n / 2; ++k)
      {
        // Binomial coefficient C(n, 2k)
        double binom = 1.0;
        for (unsigned int j = 0; j < 2 * k; ++j)
          binom *= static_cast<double>(n - j) / static_cast<double>(j + 1);

        // Double factorial (2k-1)!! with (-1)!! = 1 by convention
        double dfact = 1.0;
        for (unsigned int j = 1; j < 2 * k; j += 2)
          dfact *= j;

        moment += binom * dfact * std::pow(sigma, 2.0 * k) *
                  std::pow(mu, static_cast<double>(n - 2 * k));
      }

    return moment * std::exp((static_cast<double>(n) + alpha) * t);
  }

  /// Mean of the Gaussian distribution (phase space)
  const double mu;
  /// Standard deviation of the Gaussian distribution (phase space)
  const double sigma;
  /// Multiplication rate in the source term
  const double alpha;
};