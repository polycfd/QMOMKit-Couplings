#include "reference.hpp"

#include "exact_solution.hpp"

#include <algorithm>

namespace hyqmom_transport
{

GroupReference::GroupReference(const Problem& problem, const Mesh& mesh) : mesh_(mesh)
{
    const bool single_velocities =
        std::all_of(problem.populations.begin(), problem.populations.end(), [](const Population& p) { return p.width == 0.0; });
    if (!single_velocities)
    {
        return;
    }
    // The number density of each group at t = 0 is that of the exact solution of the group alone.
    for (const Population& population : problem.populations)
    {
        Problem single = problem;
        single.populations = {population};
        const std::vector<std::vector<double>> initial = exact_moments(single, mesh, 0.0, 1);
        std::vector<double> density(mesh.cells);
        for (std::size_t i = 0; i < mesh.cells; ++i)
        {
            density[i] = initial[i][0];
        }
        velocities_.push_back(population.velocity);
        densities_.push_back(density);
    }
}

void GroupReference::advance(const double dt)
{
    for (std::size_t p = 0; p < densities_.size(); ++p)
    {
        const double velocity = velocities_[p];
        std::vector<double>& density = densities_[p];
        std::vector<double> flux(mesh_.cells + 1);
        for (std::size_t f = 0; f <= mesh_.cells; ++f)
        {
            flux[f] = velocity * (velocity > 0.0 ? density[mesh_.left_of(f)] : density[mesh_.right_of(f)]);
        }
        for (std::size_t i = 0; i < mesh_.cells; ++i)
        {
            density[i] -= dt / mesh_.dx * (flux[i + 1] - flux[i]);
        }
    }
}

std::vector<std::vector<double>> GroupReference::moments(const std::size_t count) const
{
    if (densities_.empty())
    {
        return {};
    }
    std::vector<std::vector<double>> moments(mesh_.cells, std::vector<double>(count, 0.0));
    for (std::size_t p = 0; p < densities_.size(); ++p)
    {
        for (std::size_t i = 0; i < mesh_.cells; ++i)
        {
            double power = densities_[p][i];
            for (std::size_t k = 0; k < count; ++k)
            {
                moments[i][k] += power;
                power *= velocities_[p];
            }
        }
    }
    return moments;
}

} // namespace hyqmom_transport
