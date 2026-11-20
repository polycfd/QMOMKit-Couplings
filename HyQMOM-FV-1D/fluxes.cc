#include "fluxes.hpp"

#include <algorithm>
#include <limits>

namespace hyqmom_transport
{

namespace
{

// Adds to `flux` the moments that the weighted velocities of `cell` carry
// through a face: those moving right when `rightward` is true, and those
// moving left otherwise.
void add_kinetic_flux(const CellFluxData& cell, const bool rightward, std::vector<double>& flux)
{
    for (std::size_t j = 0; j < cell.velocities.size(); ++j)
    {
        const double velocity = cell.velocities[j];
        if ((rightward && velocity > 0.0) || (!rightward && velocity < 0.0))
        {
            double term = cell.weights[j] * velocity;
            for (double& component : flux)
            {
                component += term;
                term *= velocity;
            }
        }
    }
}

} // namespace

std::vector<double> kinetic_flux(const CellFluxData& left, const CellFluxData& right, const std::size_t count)
{
    std::vector<double> flux(count, 0.0);
    add_kinetic_flux(left, true, flux);
    add_kinetic_flux(right, false, flux);
    return flux;
}

WaveSpeeds local_wave_speeds(const CellFluxData& left, const CellFluxData& right)
{
    WaveSpeeds speeds{std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    for (const CellFluxData* cell : {&left, &right})
    {
        if (!cell->empty)
        {
            speeds.slowest = std::min(speeds.slowest, cell->slowest);
            speeds.fastest = std::max(speeds.fastest, cell->fastest);
        }
    }
    return speeds;
}

std::vector<double> hll_flux(const CellFluxData& left, const CellFluxData& right, const std::vector<double>& left_moments,
                             const std::vector<double>& right_moments, const WaveSpeeds speeds)
{
    const std::size_t count = left_moments.size();
    const std::vector<double> zero(count, 0.0);
    if (left.empty && right.empty)
    {
        return zero;
    }
    const std::vector<double>& left_flux = left.empty ? zero : left.physical_flux;
    const std::vector<double>& right_flux = right.empty ? zero : right.physical_flux;
    if (speeds.slowest >= 0.0)
    {
        return left_flux;
    }
    if (speeds.fastest <= 0.0)
    {
        return right_flux;
    }
    std::vector<double> flux(count);
    for (std::size_t k = 0; k < count; ++k)
    {
        flux[k] = (speeds.fastest * left_flux[k] - speeds.slowest * right_flux[k] +
                   speeds.slowest * speeds.fastest * (right_moments[k] - left_moments[k])) /
                  (speeds.fastest - speeds.slowest);
    }
    return flux;
}

} // namespace hyqmom_transport
