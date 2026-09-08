// Physical-space half of the Lie splitting: transport every moment with the
// same flux,
//
//     dM_k/dt + div( u M_k ) = 0
//
// The moments are transported independently, but they must not be discretised
// independently: a face value that is a convex combination of cell values
// keeps the moment *vector* inside the realizable cone, because that cone is
// convex, whereas a per-moment flux limiter distorts the ratios between the
// moments and can break realizability even when every moment is individually
// in range. This is the finite-volume counterpart of the single shared
// Zhang-Shu factor used by the deal.II program. Keeping div(phi,M) on `upwind`
// in fvSchemes is therefore the realizability-safe choice; a higher-order
// limited scheme is more accurate but has to be watched.

{
    Info<< "  Physical space: advecting " << nMoments << " moments" << endl;

    forAll(M, k)
    {
        fvScalarMatrix MEqn
        (
            fvm::ddt(M[k])
          + fvm::div(phi, M[k])
        );

        MEqn.solve();
    }
}
