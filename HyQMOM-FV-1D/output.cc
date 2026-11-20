#include "output.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <locale>
#include <sstream>

namespace hyqmom_transport
{

namespace
{

std::string time_label(const double time)
{
    std::ostringstream label;
    label << std::fixed << std::setprecision(3) << time;
    return label.str();
}

// The relative error of moment k over the cells, in the L1 norm (power 1) or the L2 norm (power 2), normalized by the same norm of `exact`.
double relative_error(const std::vector<std::vector<double>>& moments, const std::vector<std::vector<double>>& exact, const std::size_t k,
                      const int power)
{
    double difference = 0.0;
    double norm = 0.0;
    for (std::size_t i = 0; i < moments.size(); ++i)
    {
        difference += std::pow(std::abs(moments[i][k] - exact[i][k]), power);
        norm += std::pow(std::abs(exact[i][k]), power);
    }
    return std::pow(difference / norm, 1.0 / power);
}

} // namespace

Output::Output(const Settings& settings, const Problem& problem, const Mesh& mesh, const std::vector<std::vector<double>>& initial_moments)
    : directory_(settings.output / settings.case_name), label_("N" + std::to_string(settings.node_count) + "_" + flux_label(settings)), mesh_(mesh),
      initial_m0_(summed_m0(initial_moments))
{
    std::filesystem::create_directories(directory_);
    summary_.open(directory_ / (label_ + "_summary.txt"));
    summary_.imbue(std::locale::classic());
    summary_ << std::setprecision(6);

    std::ostringstream line;
    line << "case " << settings.case_name << ", N = " << settings.node_count << ", flux " << flux_label(settings) << ", " << mesh.cells
         << " cells, dx = " << mesh.dx;
    if (problem.time_step > 0.0)
    {
        line << ", dt = " << problem.time_step;
    }
    else
    {
        line << ", Courant number " << problem.courant;
    }
    report(line.str());
}

void Output::write_profiles(const double time, const std::vector<std::vector<double>>& moments, const std::vector<std::vector<double>>& exact,
                            const std::vector<std::vector<double>>& reference) const
{
    std::ofstream file(directory_ / (label_ + "_t" + time_label(time) + ".txt"));
    file.imbue(std::locale::classic());
    file << std::setprecision(17);
    const std::size_t count = moments.front().size();
    file << "# x";
    for (const std::string& prefix : {"m", "exact_m", "reference_m"})
    {
        if (prefix == "reference_m" && reference.empty())
        {
            continue;
        }
        for (std::size_t k = 0; k < count; ++k)
        {
            file << ' ' << prefix << k;
        }
    }
    file << '\n';
    for (std::size_t i = 0; i < mesh_.cells; ++i)
    {
        file << mesh_.centers[i];
        for (const double value : moments[i])
        {
            file << ' ' << value;
        }
        for (const double value : exact[i])
        {
            file << ' ' << value;
        }
        if (!reference.empty())
        {
            for (const double value : reference[i])
            {
                file << ' ' << value;
            }
        }
        file << '\n';
    }
}

void Output::report_output_time(const double time, const std::vector<std::vector<double>>& moments, const std::vector<std::vector<double>>& exact,
                                const std::vector<std::vector<double>>& reference, const std::size_t occupied_cells,
                                const std::size_t cells_with_fewer_nodes)
{
    std::ostringstream line;
    line << "t = " << time_label(time) << ": relative L1 error m0 " << relative_error(moments, exact, 0, 1) << ", m1 "
         << relative_error(moments, exact, 1, 1);
    if (!reference.empty())
    {
        line << "; relative L1 difference to the reference m0 " << relative_error(moments, reference, 0, 1) << ", m1 "
             << relative_error(moments, reference, 1, 1);
    }
    line << "; " << cells_with_fewer_nodes << " of " << occupied_cells << " occupied cells with fewer than N + 1 nodes";
    if (mesh_.periodic)
    {
        line << "; relative change of the total m0 " << (summed_m0(moments) - initial_m0_) / initial_m0_;
    }
    report(line.str());
}

void Output::report_final(const double time, const std::vector<std::vector<double>>& moments, const std::vector<std::vector<double>>& exact,
                          const std::size_t steps, const double largest_speed, const double largest_courant, const double seconds)
{
    std::ostringstream errors;
    errors << "relative L2 error at t = " << time_label(time) << ":";
    for (std::size_t k = 0; k < moments.front().size(); ++k)
    {
        errors << " m" << k << " " << relative_error(moments, exact, k, 2);
    }
    report(errors.str());

    std::ostringstream cost;
    cost << steps << " steps, largest speed " << largest_speed << ", largest Courant number from the speeds " << largest_courant << ", " << seconds
         << " s";
    report(cost.str());
}

void Output::report(const std::string& line)
{
    std::cout << line << '\n';
    summary_ << line << '\n';
}

double Output::summed_m0(const std::vector<std::vector<double>>& moments) const
{
    double sum = 0.0;
    for (const std::vector<double>& cell_moments : moments)
    {
        sum += cell_moments[0] * mesh_.dx;
    }
    return sum;
}

} // namespace hyqmom_transport
