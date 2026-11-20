// hyqmomTransport1D
//
// Free transport of the velocity moments of a number density in one space
// dimension, with the fluxes in physical space built from QMOMKit's HyQMOM.
//
// The number density f(x, xi, t) of particles at position x with velocity xi
// obeys the kinetic equation without forces,
//
//     df/dt + xi df/dx = 0,
//
// whose solution moves every particle with its own velocity:
// f(x, xi, t) = f(x - xi t, xi, 0). The code transports the raw moments
// m_k = int xi^k f dxi, k = 0, ..., 2N, of the number density in each cell,
//
//     dm_k/dt + dm_{k+1}/dx = 0,
//
// with first-order finite volumes and forward Euler. The flux of m_{2N} is the
// moment m_{2N+1}, which the transported moments do not contain. HyQMOM
// (Fox & Laurent 2022) supplies it: QMOMKit reconstructs N + 1 weighted
// points from m_0, ..., m_{2N}, and the moment m_{2N+1} of these weighted
// points closes the system. In every step, QMOMKit creates a state from the
// moments of each cell, the moments of the state replace those of the cell,
// and the code builds the face fluxes from that state in one of three ways:
//
//   nodes   The kinetic flux of Desjardins, Fox & Villedieu (2008, eqs. 33 to
//           36). The weighted points of the left cell that move right and
//           those of the right cell that move left carry their moments
//           through the face.
//   speeds  The same kinetic flux, with the 2N + 1 characteristic speeds and
//           their weights in place of the weighted points (Patel, Desjardins
//           & Fox 2019, eqs. 13 and 14).
//   hll     The HLL flux of Fox & Laurent (2022, Section 6). The physical flux
//           of each cell is (m_1, ..., m_{2N}, m_{2N+1}), and the wave speeds
//           come from the slowest and fastest characteristic speeds.
//
// There are no forces and no collisions, so nothing changes the velocities of
// the particles, and the code makes no call that advances the phase space.
//
// Usage:
//   hyqmomTransport1D --case <case> --nodes <N> --flux nodes|speeds|hll
//                     [--hll-speeds local|global] [--cells <n>]
//                     [--courant <c>] [--velocity-width <s>] [--output <dir>]
//
// The cases are lm-mono, lm-bi, lm-tri, lm-tri-gauss and lm-crossing-a, -b,
// -c (Li & Marchisio 2022, Sections 3.1 and 3.2), and fl-riemann (Fox &
// Laurent 2022, Section 6). The README describes them.

#include <qmomkit/qmomkit.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

constexpr double infinity = std::numeric_limits<double>::infinity();
constexpr double pi = 3.14159265358979323846;

// ---------------------------------------------------------------------------
// Settings from the command line
// ---------------------------------------------------------------------------

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
    std::size_t cells = 0;   // zero keeps the number of cells of the case
    double courant = 0.0;    // zero keeps the time step of the case
    double velocity_width = 0.02;
    std::filesystem::path output = "output";
};

Settings read_settings(const int argc, char* argv[])
{
    Settings settings;
    for (int i = 1; i < argc; i += 2)
    {
        const std::string option = argv[i];
        if (i + 1 >= argc)
        {
            throw std::invalid_argument("The option " + option + " needs a value.");
        }
        const std::string value = argv[i + 1];
        if (option == "--case")
        {
            settings.case_name = value;
        }
        else if (option == "--nodes")
        {
            settings.node_count = std::stoul(value);
        }
        else if (option == "--flux")
        {
            if (value == "nodes")
            {
                settings.flux = FluxType::nodes;
            }
            else if (value == "speeds")
            {
                settings.flux = FluxType::speeds;
            }
            else if (value == "hll")
            {
                settings.flux = FluxType::hll;
            }
            else
            {
                throw std::invalid_argument("The flux must be nodes, speeds or hll.");
            }
        }
        else if (option == "--hll-speeds")
        {
            if (value != "local" && value != "global")
            {
                throw std::invalid_argument("The HLL speeds must be local or global.");
            }
            settings.global_hll_speeds = value == "global";
        }
        else if (option == "--cells")
        {
            settings.cells = std::stoul(value);
        }
        else if (option == "--courant")
        {
            settings.courant = std::stod(value);
        }
        else if (option == "--velocity-width")
        {
            settings.velocity_width = std::stod(value);
        }
        else if (option == "--output")
        {
            settings.output = value;
        }
        else
        {
            throw std::invalid_argument("Unknown option " + option + ".");
        }
    }
    if (settings.case_name.empty() || settings.node_count == 0)
    {
        throw std::invalid_argument("Usage: hyqmomTransport1D --case <case> --nodes <N> --flux nodes|speeds|hll "
                                    "[--hll-speeds local|global] [--cells <n>] [--courant <c>] "
                                    "[--velocity-width <s>] [--output <dir>]");
    }
    return settings;
}

std::string flux_label(const Settings& settings)
{
    switch (settings.flux)
    {
    case FluxType::nodes:
        return "nodes";
    case FluxType::speeds:
        return "speeds";
    case FluxType::hll:
        return settings.global_hll_speeds ? "hll-global" : "hll-local";
    }
    return "";
}

// ---------------------------------------------------------------------------
// The cases
// ---------------------------------------------------------------------------

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

std::vector<double> evenly_spaced_times(const double step, const std::size_t count)
{
    std::vector<double> times;
    for (std::size_t i = 1; i <= count; ++i)
    {
        times.push_back(step * static_cast<double>(i));
    }
    return times;
}

Problem make_problem(const Settings& settings)
{
    Problem problem;
    const std::string& name = settings.case_name;

    // Li & Marchisio (2022), Section 3.1: particles on 0.3 <= x <= 0.5 in a
    // periodic domain [0, 1] with 600 cells and a time step of 1e-4 s, with
    // one, two or three velocities. Section 3.2 gives no numbers for its
    // crossing groups; the blocks, velocities and times below are read off
    // Figs. 5 and 6 of the paper.
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
    // Fox & Laurent (2022), Section 6: a Gaussian distribution of velocities
    // with m_0 = 1 and variance 1/3 everywhere, whose mean jumps from +1 for
    // x < 0 to -1 for x >= 0. Both states extend beyond the domain, so the
    // boundaries let the uniform states flow in and out.
    else if (name == "fl-riemann")
    {
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

// ---------------------------------------------------------------------------
// The exact solution
// ---------------------------------------------------------------------------

// The integrals int_lower^upper xi^k N(xi; mean, width^2) dxi, k = 0, ...,
// count - 1, of a normal distribution of velocities over [lower, upper].
// With z = (xi - mean) / width, the integrals P_i of z^i over the standard
// normal density phi follow from P_0, P_1 and an integration by parts:
// P_i = a^{i-1} phi(a) - b^{i-1} phi(b) + (i - 1) P_{i-2}, where a and b are
// the bounds in z. The binomial expansion of xi^k = (mean + width z)^k then
// gives the integrals in xi.
std::vector<double> normal_moments_over(const double mean, const double width, const double lower, const double upper,
                                        const std::size_t count)
{
    const double a = (lower - mean) / width;
    const double b = (upper - mean) / width;
    const auto density = [](const double z) { return std::isinf(z) ? 0.0 : std::exp(-0.5 * z * z) / std::sqrt(2.0 * pi); };
    // The upper tail Q(z) = int_z^inf phi, evaluated so that a difference of
    // two tails does not lose digits on either side of the mean.
    const auto upper_tail = [](const double z) { return 0.5 * std::erfc(z / std::sqrt(2.0)); };
    std::vector<double> standard(count, 0.0);
    if (a >= 0.0)
    {
        standard[0] = upper_tail(a) - upper_tail(b);
    }
    else if (b <= 0.0)
    {
        standard[0] = upper_tail(-b) - upper_tail(-a);
    }
    else
    {
        standard[0] = 1.0 - upper_tail(b) - upper_tail(-a);
    }
    if (count > 1)
    {
        standard[1] = density(a) - density(b);
    }
    for (std::size_t i = 2; i < count; ++i)
    {
        const double at_a = std::isinf(a) ? 0.0 : std::pow(a, static_cast<double>(i - 1)) * density(a);
        const double at_b = std::isinf(b) ? 0.0 : std::pow(b, static_cast<double>(i - 1)) * density(b);
        standard[i] = at_a - at_b + static_cast<double>(i - 1) * standard[i - 2];
    }

    std::vector<double> moments(count, 0.0);
    for (std::size_t k = 0; k < count; ++k)
    {
        double binomial = 1.0;
        for (std::size_t i = 0; i <= k; ++i)
        {
            moments[k] += binomial * std::pow(mean, static_cast<double>(k - i)) * std::pow(width, static_cast<double>(i)) * standard[i];
            binomial = binomial * static_cast<double>(k - i) / static_cast<double>(i + 1);
        }
    }
    return moments;
}

double overlap(const double a, const double b, const double c, const double d)
{
    return std::max(0.0, std::min(b, d) - std::max(a, c));
}

// The exact moments m_0, ..., m_{count-1}, averaged over the cell
// [cell_left, cell_right] at time t.
//
// Every particle keeps its velocity, so at time t the particles of a group
// with velocity xi are those that started at x - xi t. A group with a single
// velocity has simply moved by velocity * t, and its cell average follows from
// the overlap of the moved interval with the cell. For a group with normally
// distributed velocities, the particles at x are those whose velocities lie
// in [(x - right) / t, (x - left) / t]; the moments at x are the integrals of
// the normal distribution over this interval, and a composite Gauss-Legendre
// rule averages them over the cell. In a periodic domain, the copies of the
// interval shifted by whole domain lengths add their particles.
std::vector<double> exact_cell_moments(const Problem& problem, const double cell_left, const double cell_right, const double time,
                                       const std::size_t count)
{
    const double length = problem.x_max - problem.x_min;
    const double cell_width = cell_right - cell_left;
    std::vector<double> moments(count, 0.0);

    for (const Population& population : problem.populations)
    {
        // The velocity moments of one particle of the group: xi^k for a single
        // velocity, or the moments of the whole normal distribution.
        std::vector<double> velocity_moments(count, 1.0);
        if (population.width == 0.0)
        {
            for (std::size_t k = 1; k < count; ++k)
            {
                velocity_moments[k] = velocity_moments[k - 1] * population.velocity;
            }
        }
        else
        {
            velocity_moments = normal_moments_over(population.velocity, population.width, -infinity, infinity, count);
        }

        if (population.width == 0.0 || time == 0.0)
        {
            const double shift = population.width == 0.0 ? population.velocity * time : 0.0;
            const double left = population.left + shift;
            const double right = population.right + shift;
            double covered = 0.0;
            if (problem.periodic)
            {
                const long first = static_cast<long>(std::floor((cell_left - right) / length));
                const long last = static_cast<long>(std::ceil((cell_right - left) / length));
                for (long image = first; image <= last; ++image)
                {
                    covered += overlap(left + static_cast<double>(image) * length, right + static_cast<double>(image) * length, cell_left, cell_right);
                }
            }
            else
            {
                covered = overlap(left, right, cell_left, cell_right);
            }
            for (std::size_t k = 0; k < count; ++k)
            {
                moments[k] += population.density * velocity_moments[k] * covered / cell_width;
            }
            continue;
        }

        // A composite four-point Gauss-Legendre rule over eight sub-intervals.
        const double gauss_points[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
        const double gauss_weights[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};
        const std::size_t sub_intervals = 8;
        const double sub_width = cell_width / static_cast<double>(sub_intervals);
        // Beyond twelve standard deviations, the normal distribution holds no
        // particles that a double can represent.
        const double slowest = population.velocity - 12.0 * population.width;
        const double fastest = population.velocity + 12.0 * population.width;
        for (std::size_t s = 0; s < sub_intervals; ++s)
        {
            for (std::size_t g = 0; g < 4; ++g)
            {
                const double x = cell_left + sub_width * (static_cast<double>(s) + 0.5 * (1.0 + gauss_points[g]));
                const double weight = 0.5 * gauss_weights[g] * sub_width / cell_width;
                long first = 0;
                long last = 0;
                if (problem.periodic)
                {
                    first = static_cast<long>(std::floor((x - population.right - time * fastest) / length));
                    last = static_cast<long>(std::ceil((x - population.left - time * slowest) / length));
                }
                for (long image = first; image <= last; ++image)
                {
                    const double lower = std::max((x - population.right - static_cast<double>(image) * length) / time, slowest);
                    const double upper = std::min((x - population.left - static_cast<double>(image) * length) / time, fastest);
                    if (lower >= upper)
                    {
                        continue;
                    }
                    const std::vector<double> part = normal_moments_over(population.velocity, population.width, lower, upper, count);
                    for (std::size_t k = 0; k < count; ++k)
                    {
                        moments[k] += weight * population.density * part[k];
                    }
                }
            }
        }
    }
    return moments;
}

// ---------------------------------------------------------------------------
// Fluxes from the QMOMKit state of a cell
// ---------------------------------------------------------------------------

// What a cell contributes to the fluxes through its two faces, taken from the
// state that QMOMKit creates from the moments of the cell. An empty cell, whose
// m_0 is below QMOMKit's cutoff, has no state and contributes nothing.
struct CellFluxData
{
    bool empty = true;
    // The velocities and weights of the kinetic flux: the weighted points
    // reconstructed from the moments, or the characteristic speeds with their
    // weights.
    std::vector<double> velocities;
    std::vector<double> weights;
    // The physical flux (m_1, ..., m_{2N}, m_{2N+1}) of the weighted points,
    // for the HLL flux.
    std::vector<double> physical_flux;
    double slowest = 0.0;
    double fastest = 0.0;
    std::size_t active_node_count = 0;
};

std::string candidate_status_name(const qmomkit::StateCandidateStatus status)
{
    switch (status)
    {
    case qmomkit::StateCandidateStatus::success:
        return "success";
    case qmomkit::StateCandidateStatus::empty:
        return "empty";
    case qmomkit::StateCandidateStatus::invalid_input:
        return "invalid_input";
    case qmomkit::StateCandidateStatus::realizability_rejected:
        return "realizability_rejected";
    case qmomkit::StateCandidateStatus::reconstruction_failed:
        return "reconstruction_failed";
    }
    return "unknown";
}

// Creates the state of the cell's moments and replaces the cell's moments with
// those of the state. The moments of a cell and the weighted points
// reconstructed from them form one unit: once QMOMKit has reconstructed the
// weighted points, the cell's moments are the moments of these weighted
// points, `state.moment_values()`, as QMOMKit's coupling guide describes for
// moments that a host changes itself. The weighted points also carry the
// fluxes, so the transport then changes the moments of the same weighted
// points that it moves. When a cell holds only k <= N different velocities,
// its moments lie on the boundary of the realizable set, and the k weighted
// points determine m_{2k}, ..., m_{2N} from m_0, ..., m_{2k-1}. Moments kept
// apart from the weighted points would accumulate rounding errors in these
// higher moments that no flux removes.
CellFluxData cell_flux_data(const qmomkit::QmomCase& qmom_case, std::vector<double>& moments, const FluxType flux,
                            const std::size_t cell, const double x, const double time)
{
    CellFluxData data;

    // The state of the moments transported by this code, created the way a
    // host creates a state for moments it computed itself.
    const qmomkit::QmomStateCandidate candidate = qmom_case.state_candidate(qmomkit::MomentSet(moments));
    if (candidate.status == qmomkit::StateCandidateStatus::empty)
    {
        return data;
    }
    std::ostringstream where;
    where << "cell " << cell << " at x = " << x << ", t = " << time;
    if (!candidate.ok())
    {
        throw std::runtime_error("QMOMKit could not create a state in " + where.str() + " (" + candidate_status_name(candidate.status) +
                                 "): " + candidate.message);
    }
    const qmomkit::QmomState& state = *candidate.state;
    moments = state.moment_values();

    // The characteristic speeds of the state, with their weights and the
    // closure moment m_{2N+1}. Every flux needs the slowest and the fastest
    // speed, for the HLL wave speeds or for the time step.
    const qmomkit::HyqmomCharacteristicSpeeds speeds = qmomkit::compute_hyqmom_characteristic_speeds(state);
    if (!speeds.ok())
    {
        throw std::runtime_error("QMOMKit could not compute the characteristic speeds in " + where.str() + ": " + speeds.message);
    }

    data.empty = false;
    data.slowest = speeds.speeds.front();
    data.fastest = speeds.speeds.back();
    data.active_node_count = state.active_node_count();
    switch (flux)
    {
    case FluxType::nodes:
        data.velocities = state.quadrature().abscissas;
        data.weights = state.quadrature().weights;
        break;
    case FluxType::speeds:
        data.velocities = speeds.speeds;
        data.weights = speeds.weights;
        break;
    case FluxType::hll:
        // The state stores the moments of its weighted points, and the closure
        // moment is the next moment of the same weighted points.
        data.physical_flux.assign(state.moment_values().begin() + 1, state.moment_values().end());
        data.physical_flux.push_back(speeds.closure_moment);
        break;
    }
    return data;
}

// Adds to `flux` the moments that the weighted velocities of `cell` carry
// through a face: those moving right when `rightward` is true, and those
// moving left otherwise. A weight w at velocity v carries w v^{k+1} as the
// flux of m_k.
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

// The HLL flux between a left and a right cell, with the wave speeds
// `slowest` and `fastest`. The conserved variables are the transported moments
// of the cells, and the physical flux of an empty cell is zero.
std::vector<double> hll_flux(const CellFluxData& left, const CellFluxData& right, const std::vector<double>& left_moments,
                             const std::vector<double>& right_moments, const double slowest, const double fastest)
{
    const std::size_t count = left_moments.size();
    const std::vector<double> zero(count, 0.0);
    const std::vector<double>& left_flux = left.empty ? zero : left.physical_flux;
    const std::vector<double>& right_flux = right.empty ? zero : right.physical_flux;
    if (slowest >= 0.0)
    {
        return left_flux;
    }
    if (fastest <= 0.0)
    {
        return right_flux;
    }
    std::vector<double> flux(count);
    for (std::size_t k = 0; k < count; ++k)
    {
        flux[k] = (fastest * left_flux[k] - slowest * right_flux[k] + slowest * fastest * (right_moments[k] - left_moments[k])) /
                  (fastest - slowest);
    }
    return flux;
}

// ---------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------

std::string time_label(const double time)
{
    std::ostringstream label;
    label << std::fixed << std::setprecision(3) << time;
    return label.str();
}

void write_profiles(const std::filesystem::path& path, const std::vector<double>& centers, const std::vector<std::vector<double>>& moments,
                    const std::vector<std::vector<double>>& exact, const std::vector<std::vector<double>>& reference)
{
    std::ofstream file(path);
    file.imbue(std::locale::classic());
    file << std::setprecision(17);
    const std::size_t count = moments.front().size();
    file << "# x";
    for (std::size_t k = 0; k < count; ++k)
    {
        file << " m" << k;
    }
    for (std::size_t k = 0; k < count; ++k)
    {
        file << " exact_m" << k;
    }
    if (!reference.empty())
    {
        for (std::size_t k = 0; k < count; ++k)
        {
            file << " reference_m" << k;
        }
    }
    file << '\n';
    for (std::size_t i = 0; i < centers.size(); ++i)
    {
        file << centers[i];
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

// The relative error of moment k over the cells, in the L1 norm (power 1) or
// the L2 norm (power 2), normalized by the same norm of `exact`.
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

int main(int argc, char* argv[])
{
    try
    {
        const Settings settings = read_settings(argc, argv);
        const Problem problem = make_problem(settings);
        const std::size_t node_count = settings.node_count;
        const std::size_t moment_count = 2 * node_count + 1;
        const std::size_t cells = problem.cells;
        const double dx = (problem.x_max - problem.x_min) / static_cast<double>(cells);
        const bool fixed_time_step = problem.time_step > 0.0;

        std::vector<double> centers(cells);
        for (std::size_t i = 0; i < cells; ++i)
        {
            centers[i] = problem.x_min + (static_cast<double>(i) + 0.5) * dx;
        }
        const auto exact_at = [&](const double time)
        {
            std::vector<std::vector<double>> exact(cells);
            for (std::size_t i = 0; i < cells; ++i)
            {
                exact[i] = exact_cell_moments(problem, centers[i] - 0.5 * dx, centers[i] + 0.5 * dx, time, moment_count);
            }
            return exact;
        };

        // The initial moments are the exact cell averages at t = 0.
        std::vector<std::vector<double>> moments = exact_at(0.0);

        // The case sets the family, the node count and the options of the
        // states created from it. The HyQMOM family reconstructs N + 1 weighted
        // points from m_0, ..., m_{2N}. The default options are adaptive: when
        // the moments of a cell lie on the boundary of the realizable set, for
        // example because the cell holds only k <= N different velocities,
        // QMOMKit returns the k weighted points of a classical quadrature.
        qmomkit::QmomCase qmom_case;
        qmom_case.node_count = node_count;
        qmom_case.family = qmomkit::QmomFamily::hyperbolic;

        // When every group has a single velocity, each group can also be
        // transported on its own with first-order upwinding. The moments of
        // the groups added together are what the coupled run gives when its
        // weighted points are the velocities of the groups, so the difference
        // between the two separates the error of the closure from that of the
        // scheme.
        const bool has_reference =
            std::all_of(problem.populations.begin(), problem.populations.end(), [](const Population& p) { return p.width == 0.0; });
        std::vector<std::vector<double>> group_densities;
        if (has_reference)
        {
            for (const Population& population : problem.populations)
            {
                std::vector<double> density(cells);
                for (std::size_t i = 0; i < cells; ++i)
                {
                    Problem single = problem;
                    single.populations = {population};
                    density[i] = exact_cell_moments(single, centers[i] - 0.5 * dx, centers[i] + 0.5 * dx, 0.0, 1)[0];
                }
                group_densities.push_back(density);
            }
        }
        const auto reference_moments = [&]()
        {
            std::vector<std::vector<double>> reference(cells, std::vector<double>(moment_count, 0.0));
            for (std::size_t p = 0; p < problem.populations.size(); ++p)
            {
                for (std::size_t i = 0; i < cells; ++i)
                {
                    double power = group_densities[p][i];
                    for (std::size_t k = 0; k < moment_count; ++k)
                    {
                        reference[i][k] += power;
                        power *= problem.populations[p].velocity;
                    }
                }
            }
            return reference;
        };

        // The index of the cell on each side of face f, which lies between
        // cells f - 1 and f. Periodic boundaries wrap around; the other
        // boundaries copy the cell next to them.
        const auto left_of = [&](const std::size_t face) { return face == 0 ? (problem.periodic ? cells - 1 : 0) : face - 1; };
        const auto right_of = [&](const std::size_t face) { return face == cells ? (problem.periodic ? 0 : cells - 1) : face; };

        const std::filesystem::path directory = settings.output / settings.case_name;
        std::filesystem::create_directories(directory);
        const std::string label = "N" + std::to_string(node_count) + "_" + flux_label(settings);
        std::ofstream summary(directory / (label + "_summary.txt"));
        summary.imbue(std::locale::classic());
        summary << std::setprecision(6);

        const auto report = [&](const std::string& line)
        {
            std::cout << line << '\n';
            summary << line << '\n';
        };
        {
            std::ostringstream header;
            header << "case " << settings.case_name << ", N = " << node_count << ", flux " << flux_label(settings) << ", " << cells
                   << " cells, dx = " << dx;
            if (fixed_time_step)
            {
                header << ", dt = " << problem.time_step;
            }
            else
            {
                header << ", Courant number " << problem.courant;
            }
            report(header.str());
        }

        write_profiles(directory / (label + "_t" + time_label(0.0) + ".txt"), centers, moments, exact_at(0.0),
                       has_reference ? reference_moments() : std::vector<std::vector<double>>{});
        double initial_mass = 0.0;
        for (const std::vector<double>& cell_moments : moments)
        {
            initial_mass += cell_moments[0] * dx;
        }

        const auto clock_start = std::chrono::steady_clock::now();
        double time = 0.0;
        std::size_t steps = 0;
        double largest_speed = 0.0;
        double largest_courant = 0.0;
        std::vector<CellFluxData> cell_data(cells);
        std::vector<std::vector<double>> face_flux(cells + 1, std::vector<double>(moment_count, 0.0));

        for (const double output_time : problem.output_times)
        {
            while (true)
            {
                // 1. The state of every cell, created by QMOMKit from the cell's
                //    moments, whose moments replace those of the cell, and the
                //    extreme characteristic speeds.
                double slowest = infinity;
                double fastest = -infinity;
                for (std::size_t i = 0; i < cells; ++i)
                {
                    cell_data[i] = cell_flux_data(qmom_case, moments[i], settings.flux, i, centers[i], time);
                    if (!cell_data[i].empty)
                    {
                        slowest = std::min(slowest, cell_data[i].slowest);
                        fastest = std::max(fastest, cell_data[i].fastest);
                    }
                }
                const double largest_speed_now = std::max(std::abs(slowest), std::abs(fastest));
                largest_speed = std::max(largest_speed, largest_speed_now);

                // 2. The time step: fixed, or from the Courant number and the
                //    largest speed, with the last step ending on the output time.
                const double remaining = output_time - time;
                const double nominal_step = fixed_time_step ? problem.time_step : problem.courant * dx / largest_speed_now;
                if (remaining <= 1.0e-9 * nominal_step)
                {
                    break;
                }
                const double dt = remaining < nominal_step * (1.0 + 1.0e-9) ? remaining : nominal_step;
                largest_courant = std::max(largest_courant, dt * largest_speed_now / dx);

                // 3. The flux through every face.
                for (std::size_t f = 0; f <= cells; ++f)
                {
                    const std::size_t l = left_of(f);
                    const std::size_t r = right_of(f);
                    if (settings.flux == FluxType::hll)
                    {
                        double face_slowest = slowest;
                        double face_fastest = fastest;
                        if (!settings.global_hll_speeds)
                        {
                            // The extreme speeds of the two cells. An empty cell
                            // has no speeds, so only the other cell counts, and
                            // between two empty cells nothing moves.
                            face_slowest = infinity;
                            face_fastest = -infinity;
                            for (const CellFluxData* cell : {&cell_data[l], &cell_data[r]})
                            {
                                if (!cell->empty)
                                {
                                    face_slowest = std::min(face_slowest, cell->slowest);
                                    face_fastest = std::max(face_fastest, cell->fastest);
                                }
                            }
                        }
                        if (cell_data[l].empty && cell_data[r].empty)
                        {
                            std::fill(face_flux[f].begin(), face_flux[f].end(), 0.0);
                        }
                        else
                        {
                            face_flux[f] = hll_flux(cell_data[l], cell_data[r], moments[l], moments[r], face_slowest, face_fastest);
                        }
                    }
                    else
                    {
                        std::fill(face_flux[f].begin(), face_flux[f].end(), 0.0);
                        add_kinetic_flux(cell_data[l], true, face_flux[f]);
                        add_kinetic_flux(cell_data[r], false, face_flux[f]);
                    }
                }

                // 4. Forward Euler. Each face flux leaves one cell and enters the
                //    other, so the transport conserves m_0 summed over a periodic
                //    domain. Only the adaptive reconstruction, which drops weighted
                //    points of negligible weight, changes this sum, by far less
                //    than the errors of the scheme.
                for (std::size_t i = 0; i < cells; ++i)
                {
                    for (std::size_t k = 0; k < moment_count; ++k)
                    {
                        moments[i][k] -= dt / dx * (face_flux[i + 1][k] - face_flux[i][k]);
                    }
                }

                // 5. The groups of the reference, each upwinded with its own
                //    velocity and the same time step.
                for (std::size_t p = 0; p < group_densities.size(); ++p)
                {
                    const double velocity = problem.populations[p].velocity;
                    std::vector<double>& density = group_densities[p];
                    std::vector<double> flux(cells + 1);
                    for (std::size_t f = 0; f <= cells; ++f)
                    {
                        flux[f] = velocity * (velocity > 0.0 ? density[left_of(f)] : density[right_of(f)]);
                    }
                    for (std::size_t i = 0; i < cells; ++i)
                    {
                        density[i] -= dt / dx * (flux[i + 1] - flux[i]);
                    }
                }

                time += dt;
                ++steps;
            }
            time = output_time;

            // The profiles and the errors at this output time.
            const std::vector<std::vector<double>> exact = exact_at(time);
            const std::vector<std::vector<double>> reference =
                has_reference ? reference_moments() : std::vector<std::vector<double>>{};
            write_profiles(directory / (label + "_t" + time_label(time) + ".txt"), centers, moments, exact, reference);

            std::size_t occupied = 0;
            std::size_t fewer_nodes = 0;
            for (std::size_t i = 0; i < cells; ++i)
            {
                std::vector<double> cell_moments = moments[i];
                const CellFluxData data = cell_flux_data(qmom_case, cell_moments, settings.flux, i, centers[i], time);
                if (!data.empty)
                {
                    ++occupied;
                    if (data.active_node_count < node_count + 1)
                    {
                        ++fewer_nodes;
                    }
                }
            }
            double mass = 0.0;
            for (const std::vector<double>& cell_moments : moments)
            {
                mass += cell_moments[0] * dx;
            }

            std::ostringstream line;
            line << "t = " << time_label(time) << ": relative L1 error m0 " << relative_error(moments, exact, 0, 1) << ", m1 "
                 << relative_error(moments, exact, 1, 1);
            if (has_reference)
            {
                line << "; relative L1 difference to the reference m0 " << relative_error(moments, reference, 0, 1) << ", m1 "
                     << relative_error(moments, reference, 1, 1);
            }
            line << "; " << fewer_nodes << " of " << occupied << " occupied cells with fewer than N + 1 nodes";
            if (problem.periodic)
            {
                line << "; relative change of the total m0 " << (mass - initial_mass) / initial_mass;
            }
            report(line.str());
        }

        // The relative L2 error of every moment at the final time, normalized
        // as in Fox & Laurent (2022), Fig. 5.
        {
            const std::vector<std::vector<double>> exact = exact_at(time);
            std::ostringstream line;
            line << "relative L2 error at t = " << time_label(time) << ":";
            for (std::size_t k = 0; k < moment_count; ++k)
            {
                line << " m" << k << " " << relative_error(moments, exact, k, 2);
            }
            report(line.str());
        }
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - clock_start).count();
        std::ostringstream line;
        line << steps << " steps, largest speed " << largest_speed << ", largest Courant number from the speeds " << largest_courant << ", "
             << seconds << " s";
        report(line.str());
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
