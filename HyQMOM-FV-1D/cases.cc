#include "cases.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace hyqmom_transport
{

namespace
{

std::vector<double> evenly_spaced_times(const double step, const std::size_t count)
{
    std::vector<double> times;
    for (std::size_t i = 1; i <= count; ++i)
    {
        times.push_back(step * static_cast<double>(i));
    }
    return times;
}

} // namespace

Problem make_problem(const Settings& settings)
{
    Problem problem;
    const std::string& name = settings.case_name;

    // Li & Marchisio (2022), Section 3.1: particles on 0.3 <= x <= 0.5 in a periodic domain [0, 1] with 600 cells and a time step of 1e-4 s, with
    // one, two or three velocities. Section 3.2 gives no numbers for its crossing groups; the blocks, velocities and times below are read off Figs. 5
    // and 6 of the paper.
    if (name.rfind("lm-", 0) == 0)
    {
        problem.x_min = 0.0;
        problem.x_max = 1.0;
        problem.periodic = true;
        problem.cells = 600;
        problem.time_step = 1.0e-4;
        problem.output_times = {0.05, 0.5, 1.0, 1.5, 2.0, 2.5};
        if (name == "lm-mono")
        {
            problem.populations = {{0.3, 0.5, 1.0, 0.6, 0.0}};
        }
        else if (name == "lm-bi")
        {
            problem.populations = {{0.3, 0.5, 0.7, 0.6, 0.0}, {0.3, 0.5, 0.3, 0.4, 0.0}};
        }
        else if (name == "lm-tri" || name == "lm-tri-gauss")
        {
            const double width = name == "lm-tri" ? 0.0 : settings.velocity_width;
            problem.populations = {{0.3, 0.5, 0.5, 0.6, width}, {0.3, 0.5, 0.3, 0.4, width}, {0.3, 0.5, 0.2, 0.2, width}};
        }
        else if (name == "lm-crossing-a")
        {
            problem.populations = {{0.25, 0.35, 1.0, 0.5, 0.0}, {0.45, 0.55, 1.0, 0.0, 0.0}};
            problem.output_times = evenly_spaced_times(0.1, 6);
        }
        else if (name == "lm-crossing-b")
        {
            problem.populations = {{0.25, 0.35, 1.0, 0.5, 0.0}, {0.65, 0.75, 1.0, -0.5, 0.0}};
            problem.output_times = evenly_spaced_times(0.1, 8);
        }
        else if (name == "lm-crossing-c")
        {
            problem.populations = {{0.2, 0.3, 1.0, 0.5, 0.0}, {0.45, 0.55, 1.0, 0.0, 0.0}, {0.7, 0.8, 1.0, -0.5, 0.0}};
            problem.output_times = evenly_spaced_times(0.1, 10);
        }
        else
        {
            throw std::invalid_argument("Unknown case " + name + ".");
        }
    }
    // Fox & Laurent (2022), Section 6: a Gaussian distribution of velocities with m_0 = 1 and variance 1/3 everywhere, whose mean jumps from +1 for x
    // < 0 to -1 for x >= 0. Both states extend beyond the domain, so the boundaries let the uniform states flow in and out.
    else if (name == "fl-riemann")
    {
        const double infinity = std::numeric_limits<double>::infinity();
        problem.x_min = -0.5;
        problem.x_max = 0.5;
        problem.periodic = false;
        problem.cells = 2000;
        problem.courant = 0.9;
        problem.output_times = {0.05, 0.1};
        const double width = std::sqrt(1.0 / 3.0);
        problem.populations = {{-infinity, 0.0, 1.0, 1.0, width}, {0.0, infinity, 1.0, -1.0, width}};
    }
    else
    {
        throw std::invalid_argument("Unknown case " + name + ".");
    }

    if (settings.cells > 0)
    {
        problem.cells = settings.cells;
    }
    if (settings.courant > 0.0)
    {
        problem.courant = settings.courant;
        problem.time_step = 0.0;
    }
    return problem;
}

Mesh::Mesh(const Problem& problem)
    : cells(problem.cells), dx((problem.x_max - problem.x_min) / static_cast<double>(problem.cells)), periodic(problem.periodic),
      centers(problem.cells)
{
    for (std::size_t i = 0; i < cells; ++i)
    {
        centers[i] = problem.x_min + (static_cast<double>(i) + 0.5) * dx;
    }
}

std::size_t Mesh::left_of(const std::size_t face) const
{
    return face == 0 ? (periodic ? cells - 1 : 0) : face - 1;
}

std::size_t Mesh::right_of(const std::size_t face) const
{
    return face == cells ? (periodic ? 0 : cells - 1) : face;
}

} // namespace hyqmom_transport
