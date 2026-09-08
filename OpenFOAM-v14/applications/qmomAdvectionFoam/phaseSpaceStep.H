// Phase-space half of the Lie splitting. This is where QMOMKit does its work.
//
// In every cell independently the moment vector obeys the ODE system
//
//     dM_k/dt = S_k,     S_k = k * int xi^(k-1) G(xi) f dxi
//
// S_k is not computable from the moments alone. QMOMKit closes it with n Dirac
// deltas whose weights and abscissas come from the first 2n moments by the
// Wheeler algorithm, giving S_k = sum_i w_i k G(xi_i) xi_i^(k-1). The whole
// job - inversion and reconstruction of the rates - is the single call to
// evaluate_qmom_continuous_moment_rates() below.
//
// The system is integrated with the same explicit SSP-RK3 as the deal.II
// program. Note that the splitting is Lie, not Strang, so the scheme is
// first-order in time overall regardless of the order of the two halves.

if (growthActive)
{
    Info<< "  Phase space: evolving the moments with QMOMKit ("
        << nNodes << " nodes)" << endl;

    const scalar dt = runTime.deltaTValue();
    const scalar t0 = runTime.value() - dt;

    std::vector<double> m(nMoments);
    std::vector<double> m1(nMoments);
    std::vector<double> m2(nMoments);

    label nRepaired = 0;

    // A repaired cell is realizable by construction, but a later RK stage can
    // still leave the cone, so the retry has to be capped to keep the step
    // from looping forever on a hopeless cell.
    const label maxRepairAttempts = 5;

    // One QMOMKit inversion plus rate reconstruction. Returns false if the
    // moment sequence could not be inverted, leaving the caller to repair it
    // or to give up.
    const auto momentRates = [&]
    (
        const std::vector<double>& moments,
        const scalar t,
        std::vector<double>& rates,
        string& message
    ) -> bool
    {
        const qmomkit::QmomContinuousMomentRateResult result
        (
            qmomkit::evaluate_qmom_continuous_moment_rates
            (
                moments, nNodes, t, growthModel, qmomOptions
            )
        );

        if (!result.ok())
        {
            message = result.message.c_str();
            return false;
        }

        rates = result.rates;
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

    std::vector<double> rates(nMoments);
    string message;

    forAll(mesh.C(), celli)
    {
        for (label k = 0; k < nMoments; ++k)
        {
            m[k] = M[k][celli];
        }

        label repairAttempts = 0;

        // Three SSP-RK3 stages. A stage that cannot be inverted is retried
        // once on a repaired moment vector; if that is not allowed, or still
        // fails, the run stops rather than silently producing nonsense.
        for (label stage = 0; stage < 3; ++stage)
        {
            const std::vector<double>& in =
                (stage == 0) ? m : ((stage == 1) ? m1 : m2);

            const scalar tStage =
                (stage == 0) ? t0 : ((stage == 1) ? t0 + dt : t0 + 0.5*dt);

            if (!momentRates(in, tStage, rates, message))
            {
                if (repairBrokenCells && repairAttempts < maxRepairAttempts)
                {
                    // Repair the vector this stage started from and restart
                    // the whole step from it.
                    repair(m);
                    ++repairAttempts;

                    if (!momentRates(m, t0, rates, message))
                    {
                        FatalErrorInFunction
                            << "QMOMKit could not invert the moments in cell "
                            << celli << " at " << mesh.C()[celli]
                            << " even after repair: " << message
                            << exit(FatalError);
                    }

                    stage = -1;   // restart the RK3 from stage 0
                    continue;
                }

                FatalErrorInFunction
                    << "QMOMKit could not evaluate the phase-space moment "
                    << "rates in cell " << celli << " at " << mesh.C()[celli]
                    << nl << "    " << message << nl
                    << "    m_0 = " << m[0]
                    << ", repair attempts = " << repairAttempts << nl
                    << "Enable repair/repairBrokenCells or "
                    << "qmomOptions/adaptive in constant/qmomProperties to "
                    << "let the inversion degrade gracefully instead."
                    << exit(FatalError);
            }

            for (label k = 0; k < nMoments; ++k)
            {
                if (stage == 0)
                {
                    m1[k] = m[k] + dt*rates[k];
                }
                else if (stage == 1)
                {
                    m2[k] = 0.75*m[k] + 0.25*(m1[k] + dt*rates[k]);
                }
                else
                {
                    m[k] = (1.0/3.0)*m[k] + (2.0/3.0)*(m2[k] + dt*rates[k]);
                }
            }
        }

        if (repairAttempts > 0)
        {
            ++nRepaired;
        }

        for (label k = 0; k < nMoments; ++k)
        {
            M[k][celli] = m[k];
        }
    }

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
