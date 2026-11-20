#pragma once

// The settings of a run, read from the command line.

#include <cstddef>
#include <filesystem>
#include <string>

namespace hyqmom_transport
{

// The flux through the faces between cells, built from the QMOMKit state of each cell: the kinetic flux of the weighted points, the kinetic flux of
// the characteristic speeds, or the HLL flux.
enum class FluxType
{
    nodes,
    speeds,
    hll
};

struct Settings
{
    std::string case_name;
    std::size_t node_count = 0;
    FluxType flux = FluxType::nodes;
    bool global_hll_speeds = false;
    std::size_t cells = 0; // zero keeps the number of cells of the case
    double courant = 0.0;  // zero keeps the time step of the case
    double velocity_width = 0.02;
    std::filesystem::path output = "output";
};

Settings read_settings(int argc, char* argv[]);

// The name of the flux in the output: nodes, speeds, hll-local or hll-global.
std::string flux_label(const Settings& settings);

} // namespace hyqmom_transport
