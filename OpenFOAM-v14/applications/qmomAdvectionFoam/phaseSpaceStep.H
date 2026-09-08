// Phase-space half of the Lie splitting. This is where QMOMKit does its work.
//
// In every cell independently the number density is transported along the
// internal coordinate by the phase-space velocity G(xi),
//
//     df/dt + d( G(xi) f )/dxi = 0
//
// QMOMKit represents f by n Dirac deltas whose weights and abscissas come from
// the first 2n moments by the Wheeler algorithm, and advances that quadrature
// directly: the nodes are moved along their phase-space trajectories with an
// explicit SSP-RK3 and the raw moments are reconstructed from the advanced
// quadrature. The whole job is QmomCase::initial_state() to invert, then the
// single call to advance_qmom_continuous() below.
//
// Advecting the nodes rather than integrating dM_k/dt is what keeps the step
// realizable: the weights are carried unchanged and the abscissas only move,
// so the reconstructed moment vector cannot leave the realizable cone. Only
// the inversion of the incoming moments can fail, and that is what `repair`
// below is for. It also costs one inversion per cell per step instead of one
// per RK stage.
//
// Note that the splitting is Lie, not Strang, so the scheme is first-order in
// time overall regardless of the order of the two halves.

if (growthActive)
{
    Info<< "  Phase space: evolving the moments with QMOMKit ("
        << nNodes << " nodes)" << endl;

    const scalar dt = runTime.deltaTValue();
    const scalar t0 = runTime.value() - dt;

    // Node advection defaults to SSP-RK3; the inversion options come from
    // constant/qmomProperties. `scheme` is not used: with no discontinuous
    // moment source there is only the continuous stage to advance.
    qmomkit::QmomTransportOptions transportOptions;
    transportOptions.qmom_options = qmomOptions;

    std::vector<double> m(nMoments);

    label nRepaired = 0;

    // Advance one cell's moment vector over dt. Returns false, with a reason,
    // if the incoming sequence could not be inverted.
    const auto advance = [&]
    (
        std::vector<double>& moments,
        string& message
    ) -> bool
    {
        qmomkit::QmomCase qmomCase;
        qmomCase.moments = qmomkit::MomentSet(moments);
        qmomCase.node_count = static_cast<std::size_t>(nNodes);
        qmomCase.options = qmomOptions;

        try
        {
            const qmomkit::QmomTransportResult result
            (
                qmomkit::advance_qmom_continuous
                (
                    qmomCase.initial_state(),
                    t0,
                    dt,
                    growthModel,
                    transportOptions
                )
            );

            if (!result.ok())
            {
                message = result.message.c_str();
                return false;
            }

            // In adaptive mode the accepted state may carry fewer nodes than
            // asked for, and then reconstructs fewer moments than the solver
            // transports. Nothing sensible can be done with a short vector.
            if (label(result.state.moments().size()) != nMoments)
            {
                message =
                    "the accepted state holds "
                  + Foam::name(result.state.moments().size())
                  + " moments instead of " + Foam::name(nMoments)
                  + " (the adaptive inversion dropped nodes)";
                return false;
            }

            // A degraded inversion can still hand back a non-finite vector;
            // treat that like any other failure rather than storing it.
            for (label k = 0; k < nMoments; ++k)
            {
                if (!std::isfinite(result.state.moments()[k]))
                {
                    message = "the advanced moments are not finite";
                    return false;
                }
            }

            moments = result.state.moments();
        }
        catch (const std::exception& e)
        {
            // QmomCase::initial_state() reports a non-invertible sequence by
            // throwing rather than by a status flag.
            message = e.what();
            return false;
        }

        return true;
    };

    // Rebuild a cell that cannot be inverted from the nominal phase-space
    // shape, rescaled to the m_0 it already holds. Realizable by construction.
    const auto repair = [&](std::vector<double>& moments)
    {
        scalar m0 = max(moments[0], m0Floor);

        for (label k = 0; k < nMoments; ++k)
        {
            moments[k] = m0*nominalShape[k];
        }
    };

    string message;

    // QMOMKit follows default IEEE semantics: on a badly conditioned moment
    // sequence the Wheeler recurrence produces inf or nan, which the adaptive
    // inversion then detects, and degrades to fewer nodes. OpenFOAM traps those
    // operations instead, which turns a case QMOMKit handles into a SIGFPE. The
    // traps are therefore lifted for the duration of the step and the advanced
    // moments are checked for finiteness instead.
    const int trappedExceptions = fegetexcept();
    fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

    forAll(mesh.C(), celli)
    {
        for (label k = 0; k < nMoments; ++k)
        {
            m[k] = M[k][celli];
        }

        // A cell that cannot be inverted is retried once on a repaired moment
        // vector; if that is not allowed, or still fails, the run stops rather
        // than silently producing nonsense.
        if (!advance(m, message))
        {
            if (repairBrokenCells)
            {
                repair(m);
                ++nRepaired;

                if (!advance(m, message))
                {
                    FatalErrorInFunction
                        << "QMOMKit could not advance the moments in cell "
                        << celli << " at " << mesh.C()[celli]
                        << " even after repair: " << message
                        << exit(FatalError);
                }
            }
            else
            {
                FatalErrorInFunction
                    << "QMOMKit could not advance the phase-space step in "
                    << "cell " << celli << " at " << mesh.C()[celli] << nl
                    << "    " << message << nl
                    << "    m_0 = " << m[0] << nl
                    << "Enable repair/repairBrokenCells or "
                    << "qmomOptions/adaptive in constant/qmomProperties to "
                    << "let the inversion degrade gracefully instead."
                    << exit(FatalError);
            }
        }

        for (label k = 0; k < nMoments; ++k)
        {
            M[k][celli] = m[k];
        }
    }

    feclearexcept(FE_ALL_EXCEPT);
    feenableexcept(trappedExceptions);

    forAll(M, k)
    {
        M[k].correctBoundaryConditions();
    }

    reduce(nRepaired, sumOp<label>());

    if (nRepaired > 0)
    {
        Info<< "    repair fallback active in " << nRepaired << " cells"
            << endl;
    }

    Info<< "    min(M0) = " << gMin(M[0].primitiveField())
        << ", max(M0) = " << gMax(M[0].primitiveField()) << endl;
}
