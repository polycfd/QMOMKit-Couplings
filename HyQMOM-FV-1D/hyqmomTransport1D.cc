// hyqmomTransport1D
//
// Free transport of the velocity moments m_0, ..., m_{2N} of a number density in one space dimension, with first-order finite volumes and forward
// Euler, and with the fluxes built from QMOMKit's HyQMOM (Fox & Laurent 2022):
//
//     dm_k/dt + dm_{k+1}/dx = 0,    k = 0, ..., 2N.
//
// The flux of m_{2N} is m_{2N+1}, which the transported moments do not contain. In every step, QMOMKit reconstructs N + 1 weighted points from the
// moments of each cell, and returns the characteristic speeds of the closed system with their weights and the moment m_{2N+1} of the weighted points,
// which closes the system. The fluxes through the faces are built from these.
//
// This file holds the time loop and every call to QMOMKit, so that it shows on its own how a host uses HyQMOM. Everything from QMOMKit carries the
// prefix qmomkit::. The helpers of this program are in hyqmom_transport:
//
//   settings.hpp        the command line
//   cases.hpp           the cases and the mesh
//   exact_solution.hpp  the exact cell averages of the moments
//   fluxes.hpp          the kinetic and the HLL flux from what QMOMKit returns
//   reference.hpp       the groups of particles transported separately
//   output.hpp          the profiles and the summary
//
// The README describes the cases, the fluxes and the results.

#include "cases.hpp"
#include "exact_solution.hpp"
#include "fluxes.hpp"
#include "output.hpp"
#include "reference.hpp"
#include "settings.hpp"

#include <qmomkit/qmomkit.hpp>

#include <algorithm>
#include <cstddef>
#include <exception>
#include <iostream>
#include <vector>

int main(int argc, char* argv[])
{
    using namespace hyqmom_transport;

    try
    {
        const Settings settings = read_settings(argc, argv);
        const Problem problem = make_problem(settings);
        const Mesh mesh(problem);
        const std::size_t node_count = settings.node_count;
        const std::size_t moment_count = 2 * node_count + 1;

        // The moments m_0, ..., m_{2N} of every cell, which this code stores and transports. At t = 0, they are the exact cell averages.
        std::vector<std::vector<double>> moments = exact_moments(problem, mesh, 0.0, moment_count);

        // The QMOMKit case. It sets the family, the node count and the options of every state created from it, and its own moments are not used. The
        // hyperbolic family is HyQMOM: N + 1 weighted points from m_0, ..., m_{2N}. The default options are adaptive: when the moments of a cell lie
        // on the boundary of the realizable set, for example because the cell holds only k <= N different velocities, QMOMKit returns the k weighted
        // points of a classical quadrature.
        qmomkit::QmomCase qmom_case;
        qmom_case.node_count = node_count;
        qmom_case.family = qmomkit::QmomFamily::hyperbolic;

        GroupReference reference(problem, mesh);
        Output output(settings, problem, mesh, moments);
        output.write_profiles(0.0, moments, exact_moments(problem, mesh, 0.0, moment_count), reference.moments(moment_count));

        double time = 0.0;
        std::size_t next_output = 0;
        std::size_t steps = 0;
        double largest_speed = 0.0;
        double largest_courant = 0.0;
        std::vector<std::vector<double>> exact;
        std::vector<CellFluxData> cells(mesh.cells);
        std::vector<std::vector<double>> face_flux(mesh.cells + 1);

        while (true)
        {
            // 1. QMOMKit: the state of every cell, and what the fluxes need from it. The slowest and the fastest speed in the domain stay zero when
            //    nothing moves.
            double slowest = 0.0;
            double fastest = 0.0;
            for (std::size_t i = 0; i < mesh.cells; ++i)
            {
                cells[i] = CellFluxData{};

                // The state of the moments of the cell, created the way a host creates a state for moments it changed itself. A cell whose m_0 is
                // below QMOMKit's cutoff is empty and has no state, and any other failure stops the run.
                const qmomkit::QmomStateCandidate candidate = qmom_case.state_candidate(qmomkit::MomentSet(moments[i]));
                if (candidate.status == qmomkit::StateCandidateStatus::empty)
                {
                    continue;
                }
                if (!candidate.ok())
                {
                    output.report_qmomkit_failure(time, i, candidate.message);
                    return 1;
                }
                const qmomkit::QmomState& state = *candidate.state;

                // The moments of a cell and the weighted points reconstructed from them form one unit, so the moments of the weighted points, which
                // the state stores, replace those of the cell. The weighted points carry the fluxes, so the transport then changes the moments of the
                // same weighted points that it moves. For a cell that holds only k <= N different velocities, the k weighted points reproduce m_0,
                // ..., m_{2k-1} and determine the higher moments. Moments kept apart from the weighted points would accumulate rounding errors in
                // these higher moments that no flux removes.
                moments[i] = state.moment_values();

                // The characteristic speeds in ascending order, with a weight for each speed, and the closure moment m_{2N+1}.
                const qmomkit::HyqmomCharacteristicSpeeds speeds = qmomkit::compute_hyqmom_characteristic_speeds(state);
                if (!speeds.ok())
                {
                    output.report_qmomkit_failure(time, i, speeds.message);
                    return 1;
                }

                cells[i].empty = false;
                cells[i].slowest = speeds.speeds.front();
                cells[i].fastest = speeds.speeds.back();
                cells[i].active_node_count = state.active_node_count();
                switch (settings.flux)
                {
                case FluxType::nodes:
                    // The weighted points reconstructed from the moments.
                    cells[i].velocities = state.quadrature().abscissas;
                    cells[i].weights = state.quadrature().weights;
                    break;
                case FluxType::speeds:
                    // The characteristic speeds with their weights, which have the moments m_0, ..., m_{2N+1} of the weighted points.
                    cells[i].velocities = speeds.speeds;
                    cells[i].weights = speeds.weights;
                    break;
                case FluxType::hll:
                    // The physical flux (m_1, ..., m_{2N}, m_{2N+1}) of the weighted points.
                    cells[i].physical_flux.assign(state.moment_values().begin() + 1, state.moment_values().end());
                    cells[i].physical_flux.push_back(speeds.closure_moment);
                    break;
                }

                slowest = std::min(slowest, cells[i].slowest);
                fastest = std::max(fastest, cells[i].fastest);
            }
            const double largest_speed_now = std::max(-slowest, fastest);
            largest_speed = std::max(largest_speed, largest_speed_now);

            // 2. The output at every output time that the moments have reached, and the end of the run after the last one.
            const double nominal_step = problem.time_step > 0.0 ? problem.time_step : problem.courant * mesh.dx / largest_speed_now;
            while (next_output < problem.output_times.size() && problem.output_times[next_output] - time <= 1.0e-9 * nominal_step)
            {
                time = problem.output_times[next_output];
                exact = exact_moments(problem, mesh, time, moment_count);
                const std::vector<std::vector<double>> reference_moments = reference.moments(moment_count);
                output.write_profiles(time, moments, exact, reference_moments);
                output.report_output_time(time, moments, exact, reference_moments, cells);
                ++next_output;
            }
            if (next_output == problem.output_times.size())
            {
                break;
            }

            // 3. The time step: fixed, or from the Courant number and the largest speed, with the last step ending on the next output time. The
            //    fluxes keep the moments realizable only while no characteristic speed crosses more than one cell in a step.
            const double remaining = problem.output_times[next_output] - time;
            const double dt = remaining < nominal_step * (1.0 + 1.0e-9) ? remaining : nominal_step;
            const double courant = dt * largest_speed_now / mesh.dx;
            if (courant > 1.0)
            {
                output.report_courant_failure(time, courant);
                return 1;
            }
            largest_courant = std::max(largest_courant, courant);

            // 4. The flux through every face, between the cells on its left and on its right. With --hll-speeds global, the HLL wave speeds are the
            //    slowest and fastest speed in the domain, as in Fox & Laurent (2022), and otherwise those of the two cells.
            for (std::size_t f = 0; f <= mesh.cells; ++f)
            {
                const std::size_t l = mesh.left_of(f);
                const std::size_t r = mesh.right_of(f);
                if (settings.flux == FluxType::hll)
                {
                    const WaveSpeeds wave_speeds = settings.global_hll_speeds ? WaveSpeeds{slowest, fastest} : local_wave_speeds(cells[l], cells[r]);
                    face_flux[f] = hll_flux(cells[l], cells[r], moments[l], moments[r], wave_speeds);
                }
                else
                {
                    face_flux[f] = kinetic_flux(cells[l], cells[r], moment_count);
                }
            }

            // 5. Forward Euler. Each face flux leaves one cell and enters the other, so the transport is conservative.
            for (std::size_t i = 0; i < mesh.cells; ++i)
            {
                for (std::size_t k = 0; k < moment_count; ++k)
                {
                    moments[i][k] -= dt / mesh.dx * (face_flux[i + 1][k] - face_flux[i][k]);
                }
            }
            reference.advance(dt);

            time += dt;
            ++steps;
        }

        output.report_final(time, moments, exact, steps, largest_speed, largest_courant);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
