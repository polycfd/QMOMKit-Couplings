#pragma once

// The fluxes through the faces between cells, built from what QMOMKit returns for each cell. Nothing in this file calls QMOMKit: main() fills a
// CellFluxData for each cell from the state that QMOMKit creates.

#include <cstddef>
#include <vector>

namespace hyqmom_transport
{

// What a cell contributes to the fluxes through its two faces. An empty cell, whose m_0 is below QMOMKit's cutoff, has no state and contributes
// nothing.
struct CellFluxData
{
    bool empty = true;
    // The velocities and weights of the kinetic flux: the weighted points reconstructed from the moments, or the characteristic speeds with their
    // weights.
    std::vector<double> velocities;
    std::vector<double> weights;
    // The physical flux (m_1, ..., m_{2N}, m_{2N+1}) of the weighted points, for the HLL flux.
    std::vector<double> physical_flux;
    // The slowest and the fastest characteristic speed.
    double slowest = 0.0;
    double fastest = 0.0;
    std::size_t active_node_count = 0;
};

// The kinetic flux of m_0, ..., m_{count-1} through the face between the cells `left` and `right` (Desjardins, Fox & Villedieu 2008, eqs. 33 to 36).
// The weighted velocities of the left cell that move right and those of the right cell that move left carry their moments through the face: a weight
// w at velocity v carries w v^{k+1} as the flux of m_k.
std::vector<double> kinetic_flux(const CellFluxData& left, const CellFluxData& right, std::size_t count);

// The slowest and the fastest wave of the HLL flux.
struct WaveSpeeds
{
    double slowest;
    double fastest;
};

// The wave speeds of the HLL flux from the two cells next to a face. An empty cell has no speeds, so only the other cell counts.
WaveSpeeds local_wave_speeds(const CellFluxData& left, const CellFluxData& right);

// The HLL flux between the cells `left` and `right` with the wave speeds `speeds` (Fox & Laurent 2022, Section 6). The conserved variables are the
// moments of the cells, and the physical flux of an empty cell is zero. Between two empty cells, nothing moves.
std::vector<double> hll_flux(const CellFluxData& left, const CellFluxData& right, const std::vector<double>& left_moments,
                             const std::vector<double>& right_moments, WaveSpeeds speeds);

} // namespace hyqmom_transport
