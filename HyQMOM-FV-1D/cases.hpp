#pragma once

// The cases and the mesh: the groups of particles at t = 0, the domain, the
// time step and the output times of each case, and the cells.

#include "settings.hpp"

#include <cstddef>
#include <vector>

namespace hyqmom_transport
{

// A group of particles that at t = 0 fills the interval [left, right] with the
// number density `density`. The velocities of the group are normally
// distributed around `velocity` with the standard deviation `width`. A width
// of zero means that every particle of the group has the velocity `velocity`.
struct Population
{
    double left;
    double right;
    double density;
    double velocity;
    double width;
};

struct Problem
{
    double x_min = 0.0;
    double x_max = 1.0;
    bool periodic = true;
    std::size_t cells = 0;
    double time_step = 0.0; // a fixed time step, or zero for a time step from the Courant number
    double courant = 0.0;
    std::vector<double> output_times;
    std::vector<Population> populations;
};

// The case named in the settings, with the number of cells and the time step
// overridden by the settings when given.
Problem make_problem(const Settings& settings);

// Equally sized cells over the domain of a case.
struct Mesh
{
    explicit Mesh(const Problem& problem);

    std::size_t cells;
    double dx;
    bool periodic;
    std::vector<double> centers;

    // The cell on each side of face f, which lies between cells f - 1 and f.
    // Periodic boundaries wrap around, and the other boundaries copy the cell
    // next to them.
    [[nodiscard]] std::size_t left_of(std::size_t face) const;
    [[nodiscard]] std::size_t right_of(std::size_t face) const;
};

} // namespace hyqmom_transport
