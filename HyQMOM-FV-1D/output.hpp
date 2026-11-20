#pragma once

// The output of a run: the profiles of the moments at the output times, and
// the summary lines, which the run prints and writes to a file.

#include "cases.hpp"
#include "settings.hpp"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace hyqmom_transport
{

class Output
{
public:
    // Creates the directory <output>/<case>/ and the summary file, and reports
    // the settings of the run. The initial moments give the summed m_0 that
    // later output times compare with.
    Output(const Settings& settings, const Problem& problem, const Mesh& mesh, const std::vector<std::vector<double>>& initial_moments);

    // Writes N<N>_<flux>_t<time>.txt with the columns x, the computed moments,
    // the exact moments and, when it exists, the reference.
    void write_profiles(double time, const std::vector<std::vector<double>>& moments, const std::vector<std::vector<double>>& exact,
                        const std::vector<std::vector<double>>& reference) const;

    // Reports at an output time the relative L1 errors of m_0 and m_1, the
    // relative L1 differences to the reference, the number of occupied cells
    // whose quadrature holds fewer than N + 1 weighted points, and, in a
    // periodic domain, the relative change of the summed m_0.
    void report_output_time(double time, const std::vector<std::vector<double>>& moments, const std::vector<std::vector<double>>& exact,
                            const std::vector<std::vector<double>>& reference, std::size_t occupied_cells,
                            std::size_t cells_with_fewer_nodes);

    // Reports at the final time the relative L2 error of every moment,
    // normalized as in Fox & Laurent (2022), Fig. 5, and the cost of the run.
    void report_final(double time, const std::vector<std::vector<double>>& moments, const std::vector<std::vector<double>>& exact,
                      std::size_t steps, double largest_speed, double largest_courant, double seconds);

private:
    void report(const std::string& line);
    [[nodiscard]] double summed_m0(const std::vector<std::vector<double>>& moments) const;

    std::filesystem::path directory_;
    std::string label_;
    const Mesh& mesh_;
    double initial_m0_;
    std::ofstream summary_;
};

} // namespace hyqmom_transport
