/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    qmomAdvectionFoam

Description
    Transport of the raw moments of a number density function coupled to
    QMOMKit for the phase-space evolution.

    A number density f(x, xi, t) over physical space x and an internal
    coordinate xi obeys the population balance equation

        df/dt + div_x( u(x) f ) + d/dxi( G(xi) f ) = 0

    Rather than f itself, its raw moments M_k = int xi^k f dxi are transported:

        dM_k/dt + div_x( u M_k ) = S_k,    S_k = k int xi^(k-1) G(xi) f dxi

    The source S_k is not computable from the moments alone. QMOMKit closes it
    with the Quadrature Method of Moments: n Dirac deltas whose weights and
    abscissas are recovered from the first 2n moments by the Wheeler algorithm,
    giving S_k = sum_i w_i k G(xi_i) xi_i^(k-1).

    The two operators are Lie-split: a finite-volume advection step for the
    whole moment set, then a cell-local SSP-RK3 integration of the phase-space
    source. Everything to do with QMOM - the inversion and the reconstruction
    of the rates - is the single call to QMOMKit in phaseSpaceStep.H.

    This is the OpenFOAM counterpart of the deal.II discontinuous Galerkin
    program in DG-QMOM/QMOM_kit, and reproduces its four reference problems.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Switch.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvcDiv.H"
#include "fvMatrices.H"

#include <qmomkit/qmomkit.hpp>

#include "qmomModels.H"

#include <vector>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        #include "advectionStep.H"
        #include "phaseSpaceStep.H"
        #include "momentErrors.H"

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
