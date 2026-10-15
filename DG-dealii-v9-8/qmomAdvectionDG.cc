/* ------------------------------------------------------------------------
 *
 * QMOM solver for 1D phase and 1D physical space configuration
 *
 * Description:
 * Solves for 1D QMOM in a 1D mesh
 * This is adapted from QMOM_kit and QMOMKit::qmom_1D_continuous_breakage
 * examples
 *
 * Author: Bruna Campos adapted from Bruno Blais, Fabian Denner
 * ------------------------------------------------------------------------
 *
 */

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_function.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <qmomkit/qmomkit.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "moments_analytical_functions.hpp"
#include "moments_source_functions.hpp"
#include "velocity_functions_physical_space.hpp"

namespace QMOM

{
  /// Types of problems that can be solved
  enum class Problem
  {
    // Gaussian in phase space, no advection, linear growth in phase space
    problem_A,
    // Pure advection, no operator in phase space
    problem_B,
    // Gaussian in phase space, Gaussian in physical space
    problem_C,
    // 2D problem with a rigid body rotation
    problem_D,
    // Gaussian in phase space, no advection, linear growth in phase space,
    // multiplication source term
    problem_E,
    // Gaussian in phase space, no advection, linear growth in phase space,
    // binary fragmentation source term
    problem_F,
  };

  using namespace dealii;

  /**
   * @brief Thread-local scratch data used by MeshWorker::mesh_loop.
   *
   * Each worker thread receives its own copy of this struct so that
   * FEValues and FEInterfaceValues objects (which are not thread-safe)
   * can be reinitialised independently on each cell/face. The copy
   * constructor rebuilds fresh FEValues objects from the same mapping,
   * finite element, quadrature, and update flags.
   */
  template <int dim>
  struct ScratchData
  {
    ScratchData(const Mapping<dim>        &mapping,
                const FiniteElement<dim>  &fe,
                const Quadrature<dim>     &quadrature,
                const Quadrature<dim - 1> &quadrature_face,
                const UpdateFlags update_flags = update_values |
                                                 update_gradients |
                                                 update_quadrature_points |
                                                 update_JxW_values,
                const UpdateFlags interface_update_flags =
                  update_values | update_gradients | update_quadrature_points |
                  update_JxW_values | update_normal_vectors)
      : fe_values(mapping, fe, quadrature, update_flags)
      , fe_interface_values(mapping,
                            fe,
                            quadrature_face,
                            interface_update_flags)
    {}

    /// Copy constructor: builds new FEValues objects from the same
    /// mapping, FE, quadrature, and update flags as @p scratch_data.
    ScratchData(const ScratchData<dim> &scratch_data)
      : fe_values(scratch_data.fe_values.get_mapping(),
                  scratch_data.fe_values.get_fe(),
                  scratch_data.fe_values.get_quadrature(),
                  scratch_data.fe_values.get_update_flags())
      , fe_interface_values(scratch_data.fe_interface_values.get_mapping(),
                            scratch_data.fe_interface_values.get_fe(),
                            scratch_data.fe_interface_values.get_quadrature(),
                            scratch_data.fe_interface_values.get_update_flags())
    {}

    /// FEValues object for volume (cell) integrals.
    FEValues<dim> fe_values;
    /// FEInterfaceValues object for face integrals between cells.
    FEInterfaceValues<dim> fe_interface_values;
  };

  /**
   * @brief Per-face contribution data passed from face workers to the
   * copier in MeshWorker::mesh_loop.
   */
  struct CopyDataFace
  {
    Vector<double>                       face_rhs;
    std::vector<types::global_dof_index> joint_dof_indices;
  };

  /**
   * @brief Per-cell contribution data assembled by cell/face workers
   * and consumed by the copier (local system solver) in
   * MeshWorker::mesh_loop.
   *
   * The cell worker fills @p cell_matrix (local mass matrix) and
   * @p cell_rhs (local flux / source term). The copier then inverts the
   * mass matrix and uses the result to update the global solution vector.
   */
  struct CopyData
  {
    /// Local mass matrix (couples DOFs within a single cell).
    FullMatrix<double> cell_matrix;

    /// Local right-hand side vector (flux contributions).
    Vector<double> cell_rhs;
    /// Global DOF indices corresponding to the local DOFs on this cell.
    std::vector<types::global_dof_index> local_dof_indices;
    /// Face contributions (unused in the current no-flow problem).
    std::vector<CopyDataFace> face_data;

    /// Reinitialise all local data structures for a new cell.
    template <class Iterator>
    void
    reinit(const Iterator &cell, unsigned int dofs_per_cell)
    {
      cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
      cell_rhs.reinit(dofs_per_cell);

      local_dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
    }
  };

  /**
   * @brief Print centered header title according to given column width
   */
  std::string
  centered_column_label(const std::string &label, const std::size_t width)
  {
    if (label.size() >= width)
      return label;

    const std::size_t padding      = width - label.size();
    const std::size_t left_padding = padding / 2;
    return std::string(left_padding, ' ') + label +
           std::string(padding - left_padding, ' ');
  }

  /**
   * @brief Class for the 1D QMOM solver
   */
  template <int dim>
  class QMOM1DPhaseSpace
  {
  public:
    /// Constructor
    QMOM1DPhaseSpace(const Problem      problem,
                     const unsigned int n_moments,
                     const double       time_end,
                     const double       time_step,
                     const double       mu,
                     const double       sigma);

    /// Setup mesh (physical space)
    void
    setup_mesh();

    /// Setup system (physical space)
    void
    setup_system();

    /// Initialize moments
    void
    initialize_moments();

    /// Setup qmom options
    void
    configure_qmom();

    /// Initialize output file for L2 error
    void
    initialize_output_files();

    // Return the vector of moments at time t for a given state (this is
    // equivalent to the RHS)
    std::vector<Vector<double>>
    get_rhs_per_stage(const unsigned int stage, const FEValues<dim> &fe_v);

    /// Physical space
    void
    time_evolution_physical_space(const unsigned int stage);

    /// Phase space
    void
    time_evolution_phase_space();

    /// Record the admissible range of every moment before an advection step.
    void
    compute_limiter_bounds();

    // Apply the Zhang-Shu maximum-principle-satisfying limiter to @p solution.
    // Returns the number of cells in which the realizability safeguard had to
    // backtrack, i.e. where boundedness alone was not enough.
    unsigned int
    apply_zhang_shu_limiter(Vector<double> &solution);

    /// Whether every quadrature point of a cell carries a moment vector that
    /// QMOMKit can invert once limited with the trial factor @p theta.
    bool
    moments_are_invertible(const std::vector<Vector<double>> &quadrature_values,
                           const std::vector<double>         &averages,
                           const double                       theta);

    /// Rebuild the cells of @p solution whose moments cannot be inverted.
    /// Returns how many cells had to be repaired.
    unsigned int
    repair_moments(Vector<double> &solution);

    // Calculate L2 error between the numerical and analytical solutions
    void
    calculate_L2_error(std::ofstream &output);

    /// Output the result to Paraview
    void
    output_results(const unsigned int iteration) const;

    /// Run the simulation
    void
    solve();

    // Variables for the physical space mesh
    Triangulation<dim>   triangulation;
    const MappingQ1<dim> mapping;
    const FESystem<dim>  fe;
    DoFHandler<dim>      dof_handler;

    const QGauss<dim>     quadrature;
    const QGauss<dim - 1> quadrature_face;

    /// Velocity field in the physical space
    std::shared_ptr<TensorFunction<1, dim>> velocity_function;

    /// Solution vector of moments at the DoFs of the physical space mesh
    Vector<double>              solution_t;
    Vector<double>              solution_tdt;
    std::vector<Vector<double>> solution_stages;

    // Type of problem solved
    const Problem problem;

    // Analytical function that describes the moments
    std::shared_ptr<Function<dim>> moments_function;
    /// Solution vector of moments at current time t
    std::vector<double> moments;
    /// Number of moments
    const unsigned int n_moments;
    /// Number of quadrature nodes used for moment inversion
    const unsigned int n_nodes;
    /// Iteration number
    unsigned int iteration_number;
    /// Current time
    double time;
    /// Final time
    const double time_end;
    /// Time step
    const double time_step;
    // Output frequency
    unsigned int output_frequency = 1;
    // Output stream for error storage
    std::ofstream output_error;

    // Distribution parameters from initial moments
    const double mu          = 10;   // mean
    const double sigma       = 1.0;  // standard deviation
    const double sigma_space = 0.25; // standard deviation (physical space)

    /// Whether the Zhang-Shu limiter is applied after every Runge-Kutta stage
    /// of the physical-space (advection) step. Set to false to recover the
    /// unlimited DG scheme.
    bool use_limiter = true;

    /// Whether the limiter additionally guarantees that the moment vector it
    /// leaves behind is realizable, backtracking the Zhang-Shu factor until it
    /// is. This is what keeps a long run from dying on an impossible moment
    /// inversion; it costs one realizability check per support point per
    /// stage, and clips a little more than boundedness alone requires.
    /// Meaningless when use_limiter is false.
    bool enforce_realizability = true;

    /// Last-resort repair of cells whose moment vector cannot be inverted at
    /// all - typically the "negative recurrence" failure, where m_0 is
    /// perfectly healthy but the sequence has lost realizability and Wheeler's
    /// b coefficients turn negative. Such a cell is rebuilt from the nominal
    /// phase-space shape, rescaled to the zeroth moment the cell already
    /// carries, which is realizable by construction.
    ///
    /// This conserves m_0 exactly but overwrites the higher moments, which had
    /// stopped meaning anything anyway. Unlike the limiter it cannot fail, so
    /// it is what stands between a long run and an abort. Off by default: a
    /// run that needs it is telling you the resolution or the time step is
    /// marginal, and that is worth seeing rather than silently patching.
    bool repair_broken_cells = false;

    /// Additionally clamp m_0 up to this value when repairing, for cells that
    /// have gone negative or are numerically empty. Zero leaves m_0 alone and
    /// only repairs the shape.
    double m0_floor = 0.;

    /// Nominal phase-space shape used by the repair, i.e. the raw moments of
    /// the initial distribution normalized to m_0 = 1. Built in setup_system().
    std::vector<double> nominal_shape;

    /// How many times the Zhang-Shu factor may be halved before it is simply
    /// dropped to zero, which falls back on the cell average.
    static constexpr unsigned int max_limiter_backtracks = 30;

    /// Number of support points of one component on a cell, i.e. the number of
    /// vertices carrying a moment vector.
    const unsigned int n_support_points;

    /// Admissible range of every moment during an advection step, frozen at
    /// the beginning of the time step by compute_limiter_bounds().
    std::vector<double> limiter_min;
    std::vector<double> limiter_max;

    /// Phase space velocity
    qmomkit::PhaseSpaceVelocity phase_space_velocity;
    /// Source term
    qmomkit::QmomMomentSource source;
    /// Moment transport options
    qmomkit::QmomPhaseSpaceEvolutionOptions transport_options;
    /// Inversion options
    qmomkit::QmomOptions qmom_options;
    /// State of moments
    qmomkit::QmomCase qmom_case;

    /// Output options
    qmomkit::OutputOptions output_options;
    /// Moment writer (txt output file)
    std::optional<qmomkit::MomentWriter> moment_writer;
  };

  template <int dim>
  QMOM1DPhaseSpace<dim>::QMOM1DPhaseSpace(const Problem      problem,
                                          const unsigned int n_moments,
                                          const double       time_end,
                                          const double       time_step,
                                          const double       mu,
                                          const double       sigma)
    : mapping()
    , fe(FE_DGQ<dim>(1), n_moments)
    , dof_handler(triangulation)
    , quadrature(fe.tensor_degree() + 1)
    , quadrature_face(fe.tensor_degree() + 1)
    , problem(problem)
    , n_moments(n_moments)
    , n_nodes(n_moments / 2)
    , iteration_number(0)
    , time(0)
    , time_end(time_end)
    , time_step(time_step)
    , mu(mu)
    , sigma(sigma)
    , n_support_points(fe.base_element(0).n_dofs_per_cell())
  {
    // Fixed node count by default; a problem may override it below. QMOMKit's
    // own default is adaptive.
    qmom_options.adaptive = false;

    // Gaussian in phase space, no advection, linear growth in phase space
    if (problem == Problem::problem_A)
      {
        // Function used to compute analytical moments
        moments_function =
          std::make_shared<ConstantPhysicalGaussianPhase<dim>>(n_moments,
                                                               mu,
                                                               sigma);
        // Advection term
        phase_space_velocity = linear_growth_velocity();
        // Source term
        source = zero_source;
        // Velocity in physical space
        velocity_function = std::make_shared<NoFlow<dim>>();
      }
    // Pure advection, no operator in phase space
    else if (problem == Problem::problem_B)
      {
        // Function used to compute analytical moments
        moments_function = std::make_shared<GaussianPhysicalConstantPhase<dim>>(
          n_moments, mu, sigma, sigma_space);
        // Advection term
        phase_space_velocity = no_growth_velocity();
        // Source term
        source = zero_source;
        // Velocity in physical space
        velocity_function = std::make_shared<ConstantXFlow<dim>>();
      }
    // Gaussian in phase space, Gaussian in physical space
    else if (problem == Problem::problem_C)
      {
        // Function used to compute analytical moments
        moments_function = std::make_shared<GaussianPhysicalGaussianPhase<dim>>(
          n_moments, mu, sigma, sigma_space);
        // Advection term
        phase_space_velocity = linear_growth_velocity();
        // Source term
        source = zero_source;
        // Velocity in physical space
        velocity_function = std::make_shared<ConstantXFlow<dim>>();
        // Output frequency
        output_frequency = 10;
      }
    // 2D problem with a rigid body rotation
    else if (problem == Problem::problem_D)
      {
        qmom_options.adaptive = true;
        Point<dim> center;
        if constexpr (dim == 2)
          center[1] = 0.25;
        // Advection term
        phase_space_velocity = linear_growth_velocity();
        // Source term
        source = zero_source;
        // Velocity in physical space
        velocity_function = std::make_shared<TaylorGreenVortex<dim>>();
        // Function used to compute analytical moments
        moments_function =
          std::make_shared<SwirlGaussianPhysicalGaussianPhase<dim>>(
            n_moments, mu, sigma, sigma_space, center, velocity_function);
        // Output frequency
        output_frequency = 10;
      }
    // Gaussian in phase space, no advection, linear growth in phase space,
    // multiplication source term
    else if (problem == Problem::problem_E)
      {
        moments_function =
          std::make_shared<ConstantPhysicalGaussianPhaseWithSource<dim>>(
            n_moments, mu, sigma, 1.0);
        // Advection term
        phase_space_velocity = linear_growth_velocity();
        // Source term
        source = multiplication_source;
        // Velocity in physical space
        velocity_function = std::make_shared<NoFlow<dim>>();
      }
    // Gaussian in phase space, no advection, linear growth in phase space,
    // binary fragmentation source term
    else if (problem == Problem::problem_F)
      {
        // Function used to compute analytical moments
        moments_function =
          std::make_shared<ConstantPhysicalGaussianPhase<dim>>(n_moments,
                                                               mu,
                                                               sigma);
        // Advection term
        phase_space_velocity = linear_growth_velocity();
        // Source term
        source = binary_fragmentation_source;
        // Velocity in physical space
        velocity_function = std::make_shared<NoFlow<dim>>();
      }
  }

  /**
   * @brief Setup of the physical space mesh
   */
  template <int dim>
  void
  QMOM1DPhaseSpace<dim>::setup_mesh()
  {
    // Mesh for problems A, E and F (no advection)
    if (problem == Problem::problem_A || problem == Problem::problem_E ||
        problem == Problem::problem_F)
      {
        GridGenerator::hyper_cube(triangulation, -1, 1, true);
        std::vector<GridTools::PeriodicFacePair<
          typename Triangulation<dim>::cell_iterator>>
          periodicity_vector;
        // Apply periodic boundary conditions
        GridTools::collect_periodic_faces(
          triangulation, 0, 1, 0, periodicity_vector);
        triangulation.add_periodicity(periodicity_vector);
      }
    // Mesh for problems B (pure advection) and C (advection + phase space)
    else if (problem == Problem::problem_B || problem == Problem::problem_C)
      {
        if constexpr (dim == 1)
          GridGenerator::subdivided_hyper_cube(triangulation, 50, -1, 1, true);
        if constexpr (dim == 2)
          GridGenerator::subdivided_hyper_rectangle(
            triangulation, {100, 1}, {-1, -1}, {1, 1}, true);
        std::vector<GridTools::PeriodicFacePair<
          typename Triangulation<dim>::cell_iterator>>
          periodicity_vector;
        GridTools::collect_periodic_faces(
          triangulation, 0, 1, 0, periodicity_vector);
        triangulation.add_periodicity(periodicity_vector);
      }
    // Mesh for problem D (rigid body rotation)
    else if (problem == Problem::problem_D)
      {
        if constexpr (dim == 2)
          GridGenerator::subdivided_hyper_rectangle(
            triangulation, {40, 40}, {-1, -1}, {1, 1}, false);
      }
  }

  /**
   * @brief Setup of the physical space dofhandler and initialize solution vectors
   */
  template <int dim>
  void
  QMOM1DPhaseSpace<dim>::setup_system()
  {
    // Initialize dof handler
    dof_handler.distribute_dofs(fe);

    // Initialize solution at time t+Δt and t
    solution_tdt.reinit(dof_handler.n_dofs());
    solution_t.reinit(dof_handler.n_dofs());

    // Initialize solution stages for RK 3 stages (used in physical space
    // evolution)
    solution_stages.resize(3, Vector<double>(dof_handler.n_dofs()));
  }

  /**
   * @brief Initialize the moments using a normal distribution, as well as the
   * QmomCase object (stores the initial information regarding the phase space)
   */
  template <int dim>
  void
  QMOM1DPhaseSpace<dim>::initialize_moments()
  {
    std::cout << "Initializing moments for a normal distribution" << std::endl;
    std::cout << "Parameters: mu = " << mu << ", sigma = " << sigma
              << std::endl;

    // Compute initial moments using a normal distribution
    const qmomkit::RawMomentResult shape =
      qmomkit::compute_raw_moments(qmomkit::NormalDistribution{mu, sigma},
                                   n_moments - 1,
                                   1.0);

    AssertThrow(shape.ok(),
                ExcMessage("Could not build the nominal phase-space shape: " +
                           shape.message));

    // Store initial moments
    moments = shape.moments.values();
    moments_function->set_time(0.);

    // Interpolate initial moments to solution vectors
    VectorTools::interpolate(dof_handler, *moments_function, solution_t);
    VectorTools::interpolate(dof_handler, *moments_function, solution_tdt);

    std::cout << "Initial Moments:" << std::endl;
    for (unsigned int i = 0; i < n_moments; ++i)
      std::cout << " Moment  " << i << "  :"
                << moments_function->value(Point<dim>(), i) << std::endl;

    std::cout << std::endl;

    // Initialize qmom case (gathers the moment sequence, requested node count,
    // and inversion options)
    qmom_case = qmomkit::QmomCase{qmomkit::MomentSet(moments),
                                  n_nodes,
                                  qmom_options,
                                  std::nullopt};
  }

  /**
   * @brief QMOM configuration: set opitions and moment writer
   */
  template <int dim>
  void
  QMOM1DPhaseSpace<dim>::configure_qmom()
  {
    // Define qmom output options
    output_options.directory     = "./output";
    output_options.write_moments = true;
    output_options.interval      = 0.1;

    // Initialize moment output files and write initial moments
    std::filesystem::create_directories(output_options.directory);
    moment_writer.emplace(
      qmomkit::MomentWriter::for_qmom(output_options.directory / "moments.txt",
                                      moments.size()));
    moment_writer->write(0.0, qmomkit::MomentSet(moments));
  }

  template <int dim>
  void
  QMOM1DPhaseSpace<dim>::initialize_output_files()
  {
    const std::string file_L2error = "L2_error.dat";
    output_error.open(output_options.directory / file_L2error);
    const unsigned int column_width = 14;
    output_error << centered_column_label("time", 14);

    for (unsigned int m = 0; m < n_moments; ++m)
      {
        std::ostringstream label;
        label << "M[" << m << "]";
        output_error << " " << centered_column_label(label.str(), column_width);
      }

    output_error << "\n";
  }

  /**
   * @brief Returns the RHS of the moments equation per stage. Used in the physical
   * space evolution
   *
   * @param[in] stage Current stage
   * @param[in] fe_v FE values
   * @param[out] moments_current_stage Container that will be filled with the
   * RHS values
   *
   * @return moments_current_stage Vector filled with the interpolated RHS values
   */
  template <int dim>
  std::vector<Vector<double>>
  QMOM1DPhaseSpace<dim>::get_rhs_per_stage(const unsigned int   stage,
                                           const FEValues<dim> &fe_v)
  {
    std::vector<Vector<double>> moments_current_stage(
      fe_v.n_quadrature_points, Vector<double>(n_moments));

    // TODO add equation
    if (stage == 0)
      fe_v.get_function_values(solution_t, moments_current_stage);
    // TODO add equation
    else if (stage == 1)
      fe_v.get_function_values(solution_stages[0], moments_current_stage);
    // TODO add equation
    else if (stage == 2)
      fe_v.get_function_values(solution_stages[1], moments_current_stage);
    else
      AssertThrow(false, ExcMessage("Invalid stage."));

    return moments_current_stage;
  }

  /**
   * @brief Time evolution of the moments equation in physical space
   * Here we solve only
   * d M_k/dt + \nabla_x \cdot (u M_k) = 0
   * where u is the velocity in physical space
   *
   * @param[in] stage Current stage of the RK scheme
   */
  template <int dim>
  void
  QMOM1DPhaseSpace<dim>::time_evolution_physical_space(const unsigned int stage)
  {
    using Iterator = typename DoFHandler<dim>::active_cell_iterator;

    // Work done on individual cells
    const auto cell_worker = [&](const Iterator   &cell,
                                 ScratchData<dim> &scratch_data,
                                 CopyData         &copy_data) {
      // Get variables from scratch data
      const unsigned int n_dofs =
        scratch_data.fe_values.get_fe().n_dofs_per_cell();
      copy_data.reinit(cell, n_dofs);
      scratch_data.fe_values.reinit(cell);

      const FEValues<dim>       &fe_v           = scratch_data.fe_values;
      const std::vector<double> &JxW            = fe_v.get_JxW_values();
      std::vector<Point<dim>> quadrature_points = fe_v.get_quadrature_points();

      // Initialize RHS vector at the stage
      std::vector<Vector<double>> moments_current_stage =
        get_rhs_per_stage(stage, fe_v);

      // Loop over quadrature points
      for (unsigned int q = 0; q < fe_v.n_quadrature_points; ++q)
        {
          Tensor<1, dim> velocity =
            velocity_function->value(quadrature_points[q]);

          for (unsigned int i = 0; i < n_dofs; ++i)
            {
              // We assemble a mass matrix which we will need to invert for RK
              const unsigned int component_i =
                fe.system_to_component_index(i).first;
              for (unsigned int j = 0; j < n_dofs; ++j)
                {
                  const unsigned int component_j =
                    fe.system_to_component_index(j).first;
                  // The mass matrix only applies for every moment with itself
                  if (component_i == component_j)
                    {
                      copy_data.cell_matrix(i, j) +=
                        fe_v.shape_value(i, q)   //  \phi_i
                        * fe_v.shape_value(j, q) // \phi_j
                        * JxW[q];                // dx
                    }
                }

              // We also assemble the explicit advection term in the RHS
              copy_data.cell_rhs(i) +=
                velocity * fe_v.shape_grad(i, q) * JxW[q] *
                moments_current_stage[q][component_i]; //  u . grad(\phi_i) . M
            }
        }
    };

    // Work done on individual boundaries (nothing to do because we assume
    // periodic boundaries)
    const auto boundary_worker =
      [&]([[maybe_unused]] const Iterator     &cell,
          [[maybe_unused]] const unsigned int &face_no,
          [[maybe_unused]] ScratchData<dim>   &scratch_data,
          [[maybe_unused]] CopyData           &copy_data) {};

    // Work done on internal faces
    const auto face_worker = [&](const Iterator     &cell,
                                 const unsigned int &f,
                                 const unsigned int &sf,
                                 const Iterator     &ncell,
                                 const unsigned int &nf,
                                 const unsigned int &nsf,
                                 ScratchData<dim>   &scratch_data,
                                 CopyData           &copy_data) {
      // Get variables from scratch data
      FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;
      fe_iv.reinit(cell, f, sf, ncell, nf, nsf);
      const auto        &quadrature_points = fe_iv.get_quadrature_points();
      const unsigned int n_q_points        = quadrature_points.size();

      copy_data.face_data.emplace_back();
      CopyDataFace      &copy_data_face = copy_data.face_data.back();
      const unsigned int n_dofs         = fe.n_dofs_per_cell();
      copy_data_face.face_rhs.reinit(n_dofs);

      copy_data_face.joint_dof_indices.resize(n_dofs);
      cell->get_dof_indices(copy_data_face.joint_dof_indices);

      const std::vector<double>         &JxW     = fe_iv.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();
      // Face values on the current side of the cell
      std::vector<Vector<double>> values_here(n_q_points,
                                              Vector<double>(n_moments));
      // Face values on the adjacent cell
      std::vector<Vector<double>> values_there(n_q_points,
                                               Vector<double>(n_moments));

      if (stage == 0)
        {
          fe_iv.get_fe_face_values(0).get_function_values(solution_t,
                                                          values_here);
          fe_iv.get_fe_face_values(1).get_function_values(solution_t,
                                                          values_there);
        }
      else if (stage == 1)
        {
          fe_iv.get_fe_face_values(0).get_function_values(solution_stages[0],
                                                          values_here);
          fe_iv.get_fe_face_values(1).get_function_values(solution_stages[0],
                                                          values_there);
        }
      else if (stage == 2)
        {
          fe_iv.get_fe_face_values(0).get_function_values(solution_stages[1],
                                                          values_here);
          fe_iv.get_fe_face_values(1).get_function_values(solution_stages[1],
                                                          values_there);
        }
      else
        AssertThrow(false, ExcMessage("Invalid stage."));

      // Loop over quadrature points
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          // Fluid velocity vector
          const Tensor<1, dim> velocity =
            velocity_function->value(quadrature_points[q]);
          const double velocity_normal = velocity * normals[q];

          // Get upwind velocity
          Vector<double> &values_upwind =
            velocity_normal > 0 ? values_here[q] : values_there[q];

          // Assemble boundary terms
          for (unsigned int i = 0; i < n_dofs; ++i)
            {
              const unsigned int component_i =
                fe.system_to_component_index(i).first;

              // Boundary term contribution
              copy_data_face.face_rhs(i) -=
                fe_iv.get_fe_face_values(0).shape_value(i, q) // phi_i
                * values_upwind[component_i]                  // M
                * velocity_normal                             // u \cdot n
                * JxW[q];                                     // dx
            }
        }
    };

    // Instead of having a local copy data here we have a local system solver
    // Inverts the matrix to solver the local system and update the solution
    const auto local_system_solver = [&](const CopyData &c) {
      // Allocate local matrix and vector
      FullMatrix<double> inverse_cell_matrix(c.cell_matrix);
      Vector<double>     cell_solution(c.cell_rhs);
      Vector<double>     full_rhs(c.cell_rhs.size());
      full_rhs = c.cell_rhs;

      for (auto face_data : c.face_data)
        full_rhs += face_data.face_rhs;

      // Invert the matrix
      inverse_cell_matrix.invert(c.cell_matrix);

      // Solve the system
      inverse_cell_matrix.vmult(cell_solution, full_rhs);

      for (unsigned int i = 0; i < c.local_dof_indices.size(); ++i)
        {
          if (stage == 0)
            {
              solution_stages[stage][c.local_dof_indices[i]] =
                solution_t[c.local_dof_indices[i]] +
                time_step * cell_solution[i];
            }
          else if (stage == 1)
            {
              solution_stages[stage][c.local_dof_indices[i]] =
                0.75 * solution_t[c.local_dof_indices[i]] +
                0.25 * (solution_stages[0][c.local_dof_indices[i]] +
                        time_step * cell_solution[i]);
            }
          else if (stage == 2)
            {
              solution_stages[stage][c.local_dof_indices[i]] =
                1. / 3. * solution_t[c.local_dof_indices[i]] +
                2. / 3. *
                  (solution_stages[1][c.local_dof_indices[i]] +
                   time_step * cell_solution[i]);
            }
          else
            AssertThrow(false, ExcMessage("Invalid stage."));
        }
    };


    ScratchData<dim> scratch_data(mapping, fe, quadrature, quadrature_face);
    CopyData         copy_data;

    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          local_system_solver,
                          scratch_data,
                          copy_data,
                          MeshWorker::assemble_own_cells |
                            MeshWorker::assemble_boundary_faces |
                            MeshWorker::assemble_own_interior_faces_both,
                          boundary_worker,
                          face_worker);

    // Keep the advected moments within the bounds they already occupied. This
    // is done after every stage, since the SSP-RK3 update is a convex
    // combination of forward Euler steps and each of them has to be bounded
    // for the argument to carry over to the full step.
    if (repair_broken_cells)
      {
        const unsigned int n_repaired = repair_moments(solution_stages[stage]);
        if (n_repaired > 0)
          std::cout << "   Repaired " << n_repaired
                    << " cells with non-invertible moments" << std::endl;
      }

    if (use_limiter)
      {
        const unsigned int n_backtracked =
          apply_zhang_shu_limiter(solution_stages[stage]);
        // Silent when the Zhang-Shu factor was enough, which is the normal
        // case. A steady stream of these means the moments are close to leaving
        // the realizable cone and the resolution or the time step is marginal.
        if (n_backtracked > 0)
          std::cout << "   Realizability safeguard active in " << n_backtracked
                    << " cells" << std::endl;
      }
  }

  /**
   * @brief Time evolution of the moments equation in phase space
   * Here we solve only
   * d M_k/dt = A_k
   * where A_k includes the growth term as well as other source terms
   * (nucleation, fragmentation, aggregation, etc)
   */
  template <int dim>
  void
  QMOM1DPhaseSpace<dim>::time_evolution_phase_space()
  {
    // Current step
    const double current_step = std::min(time_step, time_end - time);

    using Iterator = typename DoFHandler<dim>::active_cell_iterator;

    // Work done on individual cells
    const auto cell_worker = [&](const Iterator   &cell,
                                 ScratchData<dim> &scratch_data,
                                 CopyData         &copy_data) {
      // Get variables from scratch data
      const unsigned int n_dofs =
        scratch_data.fe_values.get_fe().n_dofs_per_cell();
      copy_data.reinit(cell, n_dofs);
      scratch_data.fe_values.reinit(cell);

      // Why create a separate fe_v instead of fe?
      const FEValues<dim>          &fe_v = scratch_data.fe_values;
      const std::vector<double>    &JxW  = fe_v.get_JxW_values();
      const std::vector<Point<dim>> quadrature_points =
        fe_v.get_quadrature_points();

      // Interpolate moments at DOFs to DG quadrature points
      std::vector<Vector<double>> local_moments(fe_v.n_quadrature_points,
                                                Vector<double>(n_moments));
      fe_v.get_function_values(solution_t, local_moments);

      std::vector<double> evolved_local_moments(fe_v.n_quadrature_points *
                                                  n_moments,
                                                0.0);

      // Loop over quadrature points
      for (unsigned int q = 0; q < fe_v.n_quadrature_points; ++q)
        {
          // Get moments at the current quadrature point
          std::vector<double> moment_slice(n_moments, 0.0);
          for (unsigned int m = 0; m < n_moments; ++m)
            moment_slice[m] = local_moments[q][m];

          // Initialize qmom state (carries the moments and reconstruction
          // policy between time steps)
          qmomkit::QmomState state =
            qmomkit::QmomCase{qmomkit::MomentSet(moment_slice),
                              n_nodes,
                              qmom_options,
                              std::nullopt}
              .initial_state();

          // Advance all moments of each quadrature point
          // (includes both continuous and discontinuous sources)
          const qmomkit::QmomPhaseSpaceEvolutionResult evolved_moments =
            qmomkit::advance_qmom_phase_space(state,
                                              time,
                                              current_step,
                                              phase_space_velocity,
                                              source,
                                              transport_options);
          if (!evolved_moments.ok())
            throw std::runtime_error(evolved_moments.message);

          // Store in [nq x n_mom]
          for (unsigned int m = 0; m < n_moments; ++m)
            evolved_local_moments[q * n_moments + m] =
              evolved_moments.state.moment_values()[m];
        }

      for (unsigned int q = 0; q < fe_v.n_quadrature_points; ++q)
        {
          // Assemble mass matrix
          for (unsigned int i = 0; i < n_dofs; ++i)
            {
              const unsigned int component_i =
                fe.system_to_component_index(i).first;
              for (unsigned int j = 0; j < n_dofs; ++j)
                {
                  const unsigned int component_j =
                    fe.system_to_component_index(j).first;

                  // This is a diagonal mass matrix (no off-diagonal terms to
                  // couple moments)
                  if (component_i == component_j)
                    {
                      // Mass matrix contribution
                      copy_data.cell_matrix(i, j) +=
                        fe_v.shape_value(i, q)   // phi_i
                        * fe_v.shape_value(j, q) // phi_j
                        * JxW[q];                // dx
                    }
                }
            }

          // Assemble RHS with evolved moments
          for (unsigned int i = 0; i < n_dofs; ++i)
            {
              const unsigned int component_i =
                fe.system_to_component_index(i).first;

              // Project the evolved local moment values back to the DG
              // coefficients before applying the inverse mass matrix.
              copy_data.cell_rhs[i] +=
                evolved_local_moments[q * n_moments + component_i] *
                fe_v.shape_value(i, q) // phi_i
                * JxW[q];              // dx
            }
        }
    };

    // Work done on individual boundaries (nothing to do because we assume
    // periodic boundaries)
    const auto boundary_worker =
      [&]([[maybe_unused]] const Iterator     &cell,
          [[maybe_unused]] const unsigned int &face_no,
          [[maybe_unused]] ScratchData<dim>   &scratch_data,
          [[maybe_unused]] CopyData           &copy_data) {};

    // No work to do on the internal faces
    const auto face_worker =
      [&]([[maybe_unused]] const Iterator     &cell,
          [[maybe_unused]] const unsigned int &f,
          [[maybe_unused]] const unsigned int &sf,
          [[maybe_unused]] const Iterator     &ncell,
          [[maybe_unused]] const unsigned int &nf,
          [[maybe_unused]] const unsigned int &nsf,
          [[maybe_unused]] ScratchData<dim>   &scratch_data,
          [[maybe_unused]] CopyData           &copy_data) {};

    // Solve system
    const auto local_system_solver = [&](const CopyData &c) {
      // Allocate local matrix and vector
      FullMatrix<double> inverse_cell_matrix(c.cell_matrix);
      Vector<double>     cell_solution(c.cell_rhs);

      // Inverse the cell matrix
      inverse_cell_matrix.invert(c.cell_matrix);

      // Solve the system
      inverse_cell_matrix.vmult(cell_solution, c.cell_rhs);

      // Store solution in the solution vectors
      for (unsigned int i = 0; i < c.local_dof_indices.size(); ++i)
        {
          solution_t[c.local_dof_indices[i]]   = cell_solution[i];
          solution_tdt[c.local_dof_indices[i]] = cell_solution[i];
        }
    };

    ScratchData<dim> scratch_data(mapping, fe, quadrature, quadrature_face);
    CopyData         copy_data;

    // Call mesh worker: loop over cells
    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          local_system_solver,
                          scratch_data,
                          copy_data,
                          MeshWorker::assemble_own_cells |
                            MeshWorker::assemble_boundary_faces |
                            MeshWorker::assemble_own_interior_faces_once,
                          boundary_worker,
                          face_worker);

    // The phase-space stages are not limited - there is nothing being
    // transported here - but their intermediate states are fed back into the
    // next stage's inversion, so they get the same repair as the advection.
    if (repair_broken_cells)
      {
        const unsigned int n_repaired = repair_moments(solution_tdt);
        if (n_repaired > 0)
          std::cout << "   repaired " << n_repaired
                    << " cells with non-invertible moments (phase space)"
                    << std::endl;
      }
  }

  /**
   * @brief Record the pointwise range of every moment before an advection step.
   *
   * FE_DGQ is a nodal element whose support points are the cell vertices, so
   * the degrees of freedom *are* the values of the solution at those points and
   * the pointwise extrema are simply the extreme degrees of freedom. Looping
   * over cells visits every degree of freedom exactly once, the element being
   * discontinuous.
   */
  template <int dim>
  void
  QMOM1DPhaseSpace<dim>::compute_limiter_bounds()
  {
    limiter_min.assign(n_moments, std::numeric_limits<double>::max());
    limiter_max.assign(n_moments, std::numeric_limits<double>::lowest());

    std::vector<types::global_dof_index> local_dof_indices(
      fe.n_dofs_per_cell());

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
          {
            const unsigned int c     = fe.system_to_component_index(i).first;
            const double       value = solution_t[local_dof_indices[i]];
            limiter_min[c]           = std::min(limiter_min[c], value);
            limiter_max[c]           = std::max(limiter_max[c], value);
          }
      }
  }

  /**
   * @brief Zhang-Shu maximum-principle-satisfying limiter.
   *
   * Following Zhang & Shu (J. Comput. Phys. 229, 2010), the polynomial of each
   * cell is squeezed towards its own cell average,
   * \f[
   *   \tilde{u}_K(x) = \bar{u}_K + \theta_K \left( u_K(x) - \bar{u}_K \right),
   *   \qquad \theta_K \in [0, 1] ,
   * \f]
   * with \f$\theta_K\f$ the largest value that keeps the polynomial inside the
   * admissible range \f$[m, M]\f$ at the critical quadrature points:
   * \f[
   *   \theta_K = \min\left( 1,
   *     \frac{M - \bar{u}_K}{\max_K u_K - \bar{u}_K},
   *     \frac{\bar{u}_K - m}{\bar{u}_K - \min_K u_K} \right) .
   * \f]
   * The rescaling leaves the cell average untouched, so it neither creates nor
   * destroys mass; it only removes the over- and undershoots that the
   * unlimited DG polynomial produces near steep gradients. Where the solution
   * is smooth and already inside the bounds, \f$\theta_K = 1\f$ and the scheme
   * is untouched, so the formal order of accuracy is preserved.
   *
   * The admissible range \f$[m, M]\f$ is the *pointwise* range of the solution
   * entering the time step, recorded by compute_limiter_bounds(). Linear
   * advection obeys a maximum principle, so it cannot take the solution outside
   * the range it already occupies, and anything beyond that range is a
   * numerical artefact. Note the bound has to be a pointwise one: using the
   * range of the cell *averages* instead would be strictly tighter than the
   * true extrema of a peaked profile and would flatten smooth maxima. Freezing
   * the range at the beginning of the time step rather than recomputing it per
   * stage avoids ratcheting the peak down stage after stage, and letting it
   * follow the solution rather than fixing it at \f$t = 0\f$ accommodates the
   * growth imposed by the phase-space step.
   *
   * The critical points at which the polynomial has to be checked are the
   * tensor-product Gauss-Lobatto points. For the FE_DGQ(1) element used here
   * those are exactly the cell vertices, and since the element is nodal its
   * degrees of freedom *are* the values at those vertices. The limiter can
   * therefore work directly on the degrees of freedom. A higher polynomial
   * degree would require evaluating the solution at a larger set of points.
   */
  template <int dim>
  unsigned int
  QMOM1DPhaseSpace<dim>::apply_zhang_shu_limiter(Vector<double> &solution)
  {
    unsigned int       n_backtracked_cells = 0;
    const unsigned int n_dofs              = fe.n_dofs_per_cell();

    FEValues<dim> fe_values(mapping,
                            fe,
                            quadrature,
                            update_values | update_JxW_values);

    std::vector<types::global_dof_index> local_dof_indices(n_dofs);
    std::vector<Vector<double>> solution_values(quadrature.size(),
                                                Vector<double>(n_moments));

    // Squeeze every cell polynomial towards its own average.
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell->get_dof_indices(local_dof_indices);

        // Cell average of every moment. It has to be the exact average for the
        // rescaling below to be conservative; QGauss(2) integrates the Q1
        // element exactly, so it is.
        fe_values.reinit(cell);
        fe_values.get_function_values(solution, solution_values);

        double              measure = 0.;
        std::vector<double> averages(n_moments, 0.);
        for (unsigned int q = 0; q < quadrature.size(); ++q)
          {
            measure += fe_values.JxW(q);
            for (unsigned int c = 0; c < n_moments; ++c)
              averages[c] += solution_values[q][c] * fe_values.JxW(q);
          }
        for (unsigned int c = 0; c < n_moments; ++c)
          averages[c] /= measure;

        // Nodal values, gathered per moment and per support point. Since the
        // element is nodal, nodal_values[c][j] is the value of moment c at
        // vertex j, and the moment vector seen at that vertex is the column j.
        std::vector<std::vector<double>> nodal_values(
          n_moments, std::vector<double>(n_support_points, 0.));
        for (unsigned int i = 0; i < n_dofs; ++i)
          {
            const std::pair<unsigned int, unsigned int> ci =
              fe.system_to_component_index(i);
            nodal_values[ci.first][ci.second] = solution[local_dof_indices[i]];
          }

        // Zhang-Shu factor. A *single* factor is used for the whole moment
        // vector, deliberately: it makes the limited state
        //   (1 - theta) * average + theta * value
        // a convex combination of two moment vectors. Limiting each moment with
        // its own factor would distort the ratios between moments and can throw
        // the vector out of the realizable cone even though every individual
        // moment stays within its bounds.
        double theta = 1.;
        for (unsigned int c = 0; c < n_moments; ++c)
          {
            const double cell_min =
              *std::min_element(nodal_values[c].begin(), nodal_values[c].end());
            const double cell_max =
              *std::max_element(nodal_values[c].begin(), nodal_values[c].end());

            if (cell_max > averages[c])
              theta = std::min(theta,
                               (limiter_max[c] - averages[c]) /
                                 (cell_max - averages[c]));
            if (cell_min < averages[c])
              theta = std::min(theta,
                               (averages[c] - limiter_min[c]) /
                                 (averages[c] - cell_min));
          }

        // Round-off can push an average a hair outside the global range.
        theta = std::clamp(theta, 0., 1.);

        // Realizability safeguard. Bounding each moment is not enough for QMOM:
        // the moment vector as a whole has to stay in the realizable cone, or
        // the inversion has no solution and the run dies. The set of realizable
        // moment vectors is convex, and the cell average is realizable whenever
        // the polynomial is, so shrinking theta towards 0 walks the state back
        // into the cone. Backtracking by halving is enough in practice, since
        // the Zhang-Shu factor already handles the common case.
        if (enforce_realizability)
          {
            unsigned int backtracks = 0;
            while (theta > 0. &&
                   !moments_are_invertible(solution_values, averages, theta))
              {
                theta =
                  (++backtracks > max_limiter_backtracks) ? 0. : 0.5 * theta;
              }
            if (backtracks > 0)
              ++n_backtracked_cells;
          }

        for (unsigned int i = 0; i < n_dofs; ++i)
          {
            const unsigned int c     = fe.system_to_component_index(i).first;
            double            &value = solution[local_dof_indices[i]];
            value = averages[c] + theta * (value - averages[c]);
          }
      }
    return n_backtracked_cells;
  }

  /**
   * @brief Can the limited moment vectors of a cell actually be inverted?
   *
   * Builds the moment vector the trial factor @p theta would leave at every
   * quadrature point, \f$(1 - \theta) \bar{u} + \theta \, u(x_q)\f$, and asks
   * QMOMKit to invert it.
   *
   * The test is deliberately the inversion itself rather than
   * qmomkit::check_realizability(). The latter answers a stricter question and
   * is the wrong gate here: for the ten-moment Gaussian sequences this code
   * transports it reports `not_realizable` even for the exact analytical
   * moments - the five-node Hankel matrix is numerically singular and its
   * kernel-consistency test trips - while compute_qmom() inverts those very
   * same moments without trouble. Gating on it would backtrack theta to zero
   * everywhere and flatten the solution onto its cell averages. What actually
   * has to hold for the run to survive is that the phase-space step can invert
   * the moments, so that is what is checked.
   *
   * The quadrature points are also exactly the points where
   * rk_step_phase_space() later performs the inversion, so the condition is
   * verified where it is needed rather than somewhere it merely implies.
   */
  template <int dim>
  bool
  QMOM1DPhaseSpace<dim>::moments_are_invertible(
    const std::vector<Vector<double>> &quadrature_values,
    const std::vector<double>         &averages,
    const double                       theta)
  {
    std::vector<double> moments(n_moments);

    for (unsigned int q = 0; q < quadrature_values.size(); ++q)
      {
        for (unsigned int c = 0; c < n_moments; ++c)
          moments[c] =
            averages[c] + theta * (quadrature_values[q][c] - averages[c]);

        if (!qmomkit::compute_qmom(qmomkit::MomentSet(moments),
                                   n_nodes,
                                   qmom_options)
               .ok())
          return false;
      }
    return true;
  }

  /**
   * @brief Rebuild the cells whose moment vector can no longer be inverted.
   *
   * This is the fallback for the "negative recurrence" failure: Wheeler builds
   * the recurrence coefficients of the orthogonal polynomials from the moment
   * sequence, and a sequence that has lost realizability produces a negative
   * \f$b_k\f$, whose square root the Jacobi matrix needs. The zeroth moment
   * is usually perfectly healthy when this happens - it is the *shape* of the
   * sequence that has been destroyed, typically by advection debris in cells
   * that carry almost nothing.
   *
   * A cell in that state is rewritten with the nominal phase-space shape,
   * rescaled to the zeroth moment the cell already holds. Scaling a realizable
   * moment vector by a positive factor keeps it realizable - the realizable set
   * is a cone - so the result can always be inverted. The zeroth moment is
   * preserved exactly, so no mass is invented; the higher moments are
   * overwritten, but they had stopped carrying information.
   *
   * Note this is a repair, not a scheme. It cannot fail, which is exactly why
   * it is the last line of defence, but a run that leans on it heavily is
   * under-resolved rather than converged.
   */
  template <int dim>
  unsigned int
  QMOM1DPhaseSpace<dim>::repair_moments(Vector<double> &solution)
  {
    unsigned int n_repaired = 0;

    const unsigned int n_dofs = fe.n_dofs_per_cell();
    FEValues<dim>      fe_values(mapping,
                                 fe,
                                 quadrature,
                                 update_values | update_JxW_values);
    std::vector<types::global_dof_index> local_dof_indices(n_dofs);
    std::vector<Vector<double>> solution_values(quadrature.size(),
                                                Vector<double>(n_moments));
    std::vector<double>         moments(n_moments);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(solution, solution_values);

        // A cell is broken as soon as one of the points the phase-space step
        // will visit cannot be inverted.
        bool broken = false;
        for (unsigned int q = 0; q < quadrature.size() && !broken; ++q)
          {
            for (unsigned int c = 0; c < n_moments; ++c)
              moments[c] = solution_values[q][c];
            broken = !qmomkit::compute_qmom(qmomkit::MomentSet(moments),
                                            n_nodes,
                                            qmom_options)
                        .ok();
          }
        if (!broken)
          continue;

        ++n_repaired;

        // Zeroth moment the cell holds, floored if asked to.
        double measure = 0.;
        double m0      = 0.;
        for (unsigned int q = 0; q < quadrature.size(); ++q)
          {
            measure += fe_values.JxW(q);
            m0 += solution_values[q][0] * fe_values.JxW(q);
          }
        m0 = std::max(m0 / measure, m0_floor);

        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < n_dofs; ++i)
          {
            const unsigned int c = fe.system_to_component_index(i).first;
            solution[local_dof_indices[i]] = m0 * nominal_shape[c];
          }
      }
    return n_repaired;
  }

  /**
   * @brief Compute and print the L2 error between the numerical and
   * analytical solutions for each moment.
   *
   * For each moment (component) @p c, the L2 error is:
   * \f[
   *   \| e_c \|_{L^2} = \sqrt{ \int_\Omega
   *     \left( M_c^h(\mathbf{x},t) - M_c(\mathbf{x},t) \right)^2 \,
   *     \mathrm{d}\mathbf{x} }
   * \f]
   * where \f$M_c^h\f$ is the DG solution and \f$M_c\f$ the analytical
   * reference.
   */
  template <int dim>
  void
  QMOM1DPhaseSpace<dim>::calculate_L2_error(std::ofstream &output)
  {
    // Set the analytical solution at the current time
    moments_function->set_time(time);

    const unsigned int precision    = 6;
    const unsigned int column_width = 14;
    output << std::scientific << std::setprecision(precision)
           << std::setw(column_width) << time;

    // Set up FEValues to evaluate the numerical solution at quadrature points
    FEValues<dim> fe_values(mapping,
                            fe,
                            quadrature,
                            update_values | update_quadrature_points |
                              update_JxW_values);

    // Accumulate the squared error integral for each moment
    std::vector<double> error_per_moment(n_moments, 0.);

    // Storage for the numerical solution evaluated at quadrature points.
    // Each entry is a Vector<double> holding the value of every component
    // (moment) at that quadrature point.
    std::vector<Vector<double>> solution_values(
      quadrature.size(), Vector<double>(fe.n_components()));

    // Loop over all cells and quadrature points to integrate (u_h - u)^2
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(solution_tdt, solution_values);
        for (unsigned int q = 0; q < quadrature.size(); ++q)
          {
            for (unsigned int c = 0; c < fe.n_components(); ++c)
              {
                const double analytical_value =
                  moments_function->value(fe_values.get_quadrature_points()[q],
                                          c);
                const double value_moment = solution_values[q][c];
                const double error =
                  Utilities::fixed_power<2>(value_moment - analytical_value);
                error_per_moment[c] += error * fe_values.JxW(q);
              }
          }
      }

    // Take the square root to obtain the L2 norm and print per moment
    std::cout << " L2 norm of the error : " << std::endl;
    std::cout << "-----------------------" << std::endl;
    for (unsigned int c = 0; c < n_moments; ++c)
      {
        error_per_moment[c] = std::sqrt(error_per_moment[c]);
        std::cout << "Moment " << c << "  : " << error_per_moment[c]
                  << std::endl;
        output << " " << std::setprecision(precision) << std::setw(column_width)
               << error_per_moment[c];
      }
    output << "\n";
    std::cout << std::endl;
  }

  template <int dim>
  void
  QMOM1DPhaseSpace<dim>::output_results(const unsigned int iteration) const
  {
    const std::string prefix = (problem == Problem::problem_A) ? "solutionA" :
                               (problem == Problem::problem_B) ? "solutionB" :
                               (problem == Problem::problem_C) ? "solutionC" :
                               (problem == Problem::problem_D) ? "solutionD" :
                               (problem == Problem::problem_E) ? "solutionE" :
                                                                 "solutionF";

    const std::string filename =
      prefix + "-" + std::to_string(iteration) + ".vtk";
    std::cout << "  Writing solution to <" << filename << '>' << std::endl;
    std::ofstream output(output_options.directory / filename);
    std::cout << std::endl;

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution_tdt, "u", DataOut<dim>::type_dof_data);

    data_out.build_patches(mapping);

    data_out.write_vtk(output);
  }

  template <int dim>
  void
  QMOM1DPhaseSpace<dim>::solve()
  {
    // Setup triangulation
    setup_mesh();

    // Setup vectors and initialize dof handler
    setup_system();

    // Initialize moments
    initialize_moments();

    // Configure qmom options
    configure_qmom();

    // Initialize file to store L2 error
    initialize_output_files();

    // Initialize time variables
    time             = 0;
    iteration_number = 0;
    output_results(iteration_number);

    // Initialize variables to track time output interval
    const double output_interval  = output_options.interval.value_or(time_end);
    double       next_output_time = output_interval;

    // Loop over time
    while (time < (time_end - 1e-12))
      {
        std::cout << "****************************" << std::endl;
        std::cout << "Time step: " << iteration_number << ", time = " << time
                  << std::endl;
        std::cout << "****************************" << std::endl;

        // Update time variables
        iteration_number++;
        time += time_step;

        if (use_limiter)
          compute_limiter_bounds();

        // Physical space evolution
        for (unsigned int stage = 0; stage < 3; ++stage)
          time_evolution_physical_space(stage);

        // Update solution vectors
        solution_tdt = solution_stages[2];
        solution_t   = solution_tdt;

        // Phase space evolution
        time_evolution_phase_space();

        // Print analytical moments (not available for problem F)
        moments_function->set_time(time);
        if (!(problem == QMOM::Problem::problem_F))
          {
            std::cout << " Analytical Moments:" << std::endl;
            std::cout << "--------------------" << std::endl;
            for (unsigned int i = 0; i < n_moments; ++i)
              std::cout << " Moment  " << i << "  :"
                        << moments_function->value(Point<dim>(), i)
                        << std::endl;

            std::cout << std::endl;

            calculate_L2_error(output_error);
          }

        const bool is_last_step = time >= (time_end - 1e-12);
        if (iteration_number % output_frequency == 0 || is_last_step)
          output_results(iteration_number);

        // Save updated moments in output file
        if (time + 1.0e-12 >= next_output_time || time >= time_end - 1.0e-12)
          {
            moment_writer->write(time, qmomkit::MomentSet(moments));

            while (next_output_time <= time + 1.0e-12)
              next_output_time += output_interval;
          }
      }
  }


} // namespace QMOM

int
main(int argc, char *argv[])
{
  try
    {
      // Select problem type
      QMOM::Problem problem;
      if (argc < 2)
        {
          std::cout << "WARNING: No problem specified. Defaulting to problem A."
                    << std::endl;
          std::cout << "Usage: " << argv[0] << " <A|B|C|D|E|F>" << std::endl;
          problem = QMOM::Problem::problem_A;
        }
      else
        {
          const std::string arg(argv[1]);
          if (arg == "A" || arg == "a")
            problem = QMOM::Problem::problem_A;
          else if (arg == "B" || arg == "b")
            problem = QMOM::Problem::problem_B;
          else if (arg == "C" || arg == "c")
            problem = QMOM::Problem::problem_C;
          else if (arg == "D" || arg == "d")
            problem = QMOM::Problem::problem_D;
          else if (arg == "E" || arg == "e")
            problem = QMOM::Problem::problem_E;
          else if (arg == "F" || arg == "f")
            problem = QMOM::Problem::problem_F;
          else
            throw std::invalid_argument("Invalid problem '" + arg +
                                        "'. Valid options are <A|B|C|D|E|F>.");

          std::cout << "---------- Problem type: " << arg << " ----------"
                    << std::endl;
        }

      /// Number of moments
      const unsigned int N = 10;
      /// Final simulation time
      const double time_end = 2;
      /// Time step
      const double time_step = 0.002;
      /// Mean used in normal distribution (generate initial moments)
      const double mu = 10;
      /// Standard deviation used in normal distribution (generate initial
      /// moments)
      const double sigma = 1;

      std::cout << "Number of moments: " << N << ", t_end = " << time_end
                << ", time_step = " << time_step << std::endl;
      std::cout << "-------------------------------------" << std::endl;

      // Solve problem
      // 2D case for problem D (TGV), 1D otherwise
      if (problem == QMOM::Problem::problem_D)
        {
          QMOM::QMOM1DPhaseSpace<2> qmom(
            problem, N, time_end, time_step, mu, sigma);
          qmom.solve();
        }
      else
        {
          QMOM::QMOM1DPhaseSpace<1> qmom(
            problem, N, time_end, time_step, mu, sigma);
          qmom.solve();
        }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}