
#include <qmomkit/qmomkit.hpp>


/**
 * @brief Zero phase-space velocity
 *
 * @param[in] t Current time (unused)
 * @param[in] xi Phase space vector (unused)
 *
 * @return Vector for zero growth velocity
 */
inline qmomkit::PhaseSpaceVelocity
no_growth_velocity()
{
  return [](const double /*t*/,
            const std::vector<double> &xi,
            std::vector<double>       &velocity,
            unsigned int) {
    for (unsigned int node = 0; node < xi.size(); ++node)
      velocity[node] = 0.0;
  };
}

/**
 * @brief Linear growth velocity
 * G(\xi) = 0.01 * \xi
 *
 * @param[in] t Current time (unused)
 * @param[in] xi Phase space vector
 *
 * @return Vector for linear growth velocity
 */
inline qmomkit::PhaseSpaceVelocity
linear_growth_velocity()
{
  return [](const double /*t*/,
            const std::vector<double> &xi,
            std::vector<double>       &velocity,
            unsigned int) {
    for (unsigned int node = 0; node < xi.size(); ++node)
      velocity[node] = 0.01 * xi[node];
  };
}

/**
 * @brief Binary fragmentation source term
 * S = rate * (2^(k-1) - 1) \xi^k, or in terms of the moments:
 * S = rate * (2^(k-1) - 1) M_k
 * There is no analytical function corresponding to this term.
 *
 * @param[in] t Current time (unused)
 * @param[in] moments Current moments stored in QMOMKit
 *
 * @return Vector for binary fragmentation source term
 */
inline qmomkit::QmomMomentRate
binary_fragmentation_source(const double /*t*/,
                            const qmomkit::QmomSourceContext &moments)
{
  // Fragmentation rate
  // TODO pass as function argument
  const double fragmentation_rate = 0.5;

  std::vector<double> source(moments.size(), 0.0);
  for (unsigned int k = 0; k < moments.size(); ++k)
    {
      const double aux = std::pow(2.0, 1.0 - k);
      source[k]        = fragmentation_rate * (aux - 1) * moments.value(k);
    }

  return qmomkit::QmomMomentRate{std::move(source)};
};

/**
 * @brief Multiplication source term
 * S = rate * \xi^k, or in terms of the moments:
 * S = rate * M_k
 * This source term still yields a simple analytical solution, implemented in
 * moments_analytical_functions.hpp
 *
 * @param[in] t Current time (unused)
 * @param[in] moments Current moments stored in QMOMKit
 *
 * @return Vector for multiplication source term
 */
inline qmomkit::QmomMomentRate
multiplication_source(const double /*t*/,
                      const qmomkit::QmomSourceContext &moments)
{
  // Multiplication rate
  // TODO pass as function argument
  const double rate = 1.0;

  std::vector<double> source(moments.size(), 0.0);
  for (unsigned int k = 0; k < moments.size(); ++k)
    source[k] = rate * moments.value(k);

  return qmomkit::QmomMomentRate{std::move(source)};
};

/**
 * @brief Zero source term
 *
 * @param[in] t Current time (unused)
 * @param[in] moments Current moments stored in QMOMKit (unused)
 *
 * @return Vector for zero source term
 */
inline qmomkit::QmomMomentRate
zero_source(const double /*t*/, const qmomkit::QmomSourceContext &moments)
{
  std::vector<double> source(moments.size(), 0.0);

  return qmomkit::QmomMomentRate{std::move(source)};
};