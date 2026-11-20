#include "exact_solution.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace hyqmom_transport
{

namespace
{

constexpr double infinity = std::numeric_limits<double>::infinity();
constexpr double pi = 3.14159265358979323846;

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
// in [(x - right) / t, (x - left) / t]. The moments at x are the integrals of
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
                    const double shift_by_image = static_cast<double>(image) * length;
                    covered += overlap(left + shift_by_image, right + shift_by_image, cell_left, cell_right);
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
                    const double shift_by_image = static_cast<double>(image) * length;
                    const double lower = std::max((x - population.right - shift_by_image) / time, slowest);
                    const double upper = std::min((x - population.left - shift_by_image) / time, fastest);
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

} // namespace

std::vector<std::vector<double>> exact_moments(const Problem& problem, const Mesh& mesh, const double time, const std::size_t count)
{
    std::vector<std::vector<double>> moments(mesh.cells);
    for (std::size_t i = 0; i < mesh.cells; ++i)
    {
        moments[i] = exact_cell_moments(problem, mesh.centers[i] - 0.5 * mesh.dx, mesh.centers[i] + 0.5 * mesh.dx, time, count);
    }
    return moments;
}

} // namespace hyqmom_transport
