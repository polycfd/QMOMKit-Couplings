#pragma once

// The exact solution of the cases: every particle keeps its velocity, so the number density at time t is the initial number density moved along the
// velocities, n(x, xi, t) = n(x - xi t, xi, 0).

#include "cases.hpp"

#include <cstddef>
#include <vector>

namespace hyqmom_transport
{

// The exact moments m_0, ..., m_{count-1} at time t, averaged over each cell of the mesh.
std::vector<std::vector<double>> exact_moments(const Problem& problem, const Mesh& mesh, double time, std::size_t count);

} // namespace hyqmom_transport
