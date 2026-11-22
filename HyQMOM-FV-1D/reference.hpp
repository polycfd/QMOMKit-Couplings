#pragma once

// A reference for the cases whose groups of particles each have a single velocity: every group is transported on its own with first-order upwinding
// and its own velocity, with the same time steps as the coupled run. The moments of the groups added together are what the coupled run gives when the
// weighted points of every cell are the velocities of the groups. The difference between the two therefore separates the error of the closure from
// that of the scheme.

#include "cases.hpp"

#include <cstddef>
#include <vector>

namespace hyqmom_transport
{

class GroupReference
{
public:
    GroupReference(const Problem& problem, const Mesh& mesh);

    // Transports every group over the time step dt.
    void advance(double dt);

    // The moments m_0, ..., m_{count-1} of the groups added together in each cell, or nothing when no reference exists.
    [[nodiscard]] std::vector<std::vector<double>> moments(std::size_t count) const;

private:
    const Mesh& mesh_;
    std::vector<double> velocities_;
    std::vector<std::vector<double>> densities_;
};

} // namespace hyqmom_transport
