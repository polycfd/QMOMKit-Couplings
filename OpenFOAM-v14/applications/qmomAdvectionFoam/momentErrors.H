// L2 error of every moment against the analytical solution, when the problem
// has one. The absolute error alone is hard to read: the high moments are
// enormous (m_9 ~ 9e16 at t = 2 for the growing Gaussian), so the relative
// error is reported next to it.

if
(
    exactSolution.active()
 && (
        runTime.timeIndex() % exactSolution.reportEvery() == 0
     || runTime.writeTime()
     || !runTime.running()          // always report the final step
    )
)
{
    const scalar t = runTime.value();

    Info<< "  L2 error against the analytical solution:" << endl;

    forAll(M, k)
    {
        scalar sumErrSqr = 0;
        scalar sumRefSqr = 0;

        forAll(mesh.C(), celli)
        {
            const scalar ref = exactSolution.value(k, mesh.C()[celli], t);
            const scalar V = mesh.V()[celli];

            sumErrSqr += V*sqr(M[k][celli] - ref);
            sumRefSqr += V*sqr(ref);
        }

        reduce(sumErrSqr, sumOp<scalar>());
        reduce(sumRefSqr, sumOp<scalar>());

        const scalar absErr = Foam::sqrt(sumErrSqr);
        const scalar relErr =
            (sumRefSqr > vSmall) ? Foam::sqrt(sumErrSqr/sumRefSqr) : 0;

        Info<< "    M" << k << ": abs = " << absErr
            << ", rel = " << relErr << endl;
    }
}
