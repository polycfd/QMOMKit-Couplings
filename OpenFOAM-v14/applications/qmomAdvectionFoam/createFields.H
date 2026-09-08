// * * * * * * * * * * * * * *  QMOM properties  * * * * * * * * * * * * * * //

Info<< "Reading qmomProperties\n" << endl;

IOdictionary qmomProperties
(
    IOobject
    (
        "qmomProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const label nMoments(qmomProperties.lookup<label>("nMoments"));

// A quadrature with n nodes is built from moments 0 ... 2n-1, so the natural
// choice is n = nMoments/2, exactly as in the deal.II program.
const label nNodes
(
    qmomProperties.lookupOrDefault<label>("nNodes", nMoments/2)
);

if (nMoments < 2 || nMoments % 2 != 0)
{
    FatalErrorInFunction
        << "nMoments must be a positive even number, not " << nMoments
        << exit(FatalError);
}

if (2*nNodes > nMoments)
{
    FatalErrorInFunction
        << "A " << nNodes << "-node quadrature needs " << 2*nNodes
        << " moments, but only " << nMoments << " are transported"
        << exit(FatalError);
}

// Numerical controls handed to every QMOMKit inversion.
//
// adaptive makes Wheeler degrade gracefully rather than fail outright: it
// drops to fewer nodes when the moment sequence is nearly degenerate, and to a
// single node once m_0 falls below rmin. That is the cure for a negative
// recurrence coefficient, which is what an unrealizable sequence produces.
// It is not free - it also fires on merely ill-conditioned sequences and
// inflates the error of the high moments - so it is off by default.
qmomkit::QmomOptions qmomOptions;
{
    const dictionary& dict = qmomProperties.subOrEmptyDict("qmomOptions");

    qmomOptions.adaptive = dict.lookupOrDefault<Switch>("adaptive", false);
    qmomOptions.rmin = dict.lookupOrDefault<scalar>("rmin", 1e-10);
    qmomOptions.eabs = dict.lookupOrDefault<scalar>("eabs", 1e-10);
    qmomOptions.cutoff = dict.lookupOrDefault<scalar>("cutoff", 1e-30);
}

// * * * * * * * * * * * * *  Phase-space growth  * * * * * * * * * * * * * * //

const dictionary& phaseSpaceDict = qmomProperties.subDict("phaseSpace");

const word growthModelType(phaseSpaceDict.lookup<word>("growthModel"));
const scalar growthRate
(
    phaseSpaceDict.lookupOrDefault<scalar>("growthRate", 1)
);

const qmomkit::PhaseSpaceVelocity growthModel
(
    qmom::newGrowthModel(growthModelType, growthRate)
);

const bool growthActive = (growthModelType != "none");

// * * * * * * * * * * *  Initial phase-space distribution  * * * * * * * * * //

const dictionary& initDict = qmomProperties.subDict("initialDistribution");

const scalar mu(initDict.lookup<scalar>("mu"));
const scalar sigma(initDict.lookup<scalar>("sigma"));
const word spatialShape(initDict.lookupOrDefault<word>("spatialShape", "uniform"));
const scalar sigmaSpace(initDict.lookupOrDefault<scalar>("sigmaSpace", 0.25));
const point blobCentre(initDict.lookupOrDefault<point>("centre", Zero));

// Raw moments E[X^k] of X ~ N(mu, sigma), from QMOMKit itself. These are both
// the initial phase-space shape and, normalised to m_0 = 1, the nominal shape
// the repair fallback rebuilds a broken cell from.
std::vector<double> phaseMoments;
{
    const qmomkit::RawMomentResult shape
    (
        qmomkit::compute_raw_moments
        (
            qmomkit::NormalDistribution{mu, sigma},
            nMoments - 1,
            1.0
        )
    );

    if (!shape.ok())
    {
        FatalErrorInFunction
            << "QMOMKit could not build the initial phase-space moments: "
            << shape.message.c_str() << exit(FatalError);
    }

    phaseMoments = shape.moments.values();
}

Info<< "Initial phase-space moments of N(" << mu << ", " << sigma << "):" << nl;
for (label k = 0; k < nMoments; ++k)
{
    Info<< "    E[X^" << k << "] = " << phaseMoments[k] << nl;
}
Info<< endl;

// * * * * * * * * * * * * * * *  Velocity field  * * * * * * * * * * * * * * //

const qmom::velocityModel velocityModel
(
    qmomProperties.subDict("velocityModel")
);

Info<< "Creating velocity field U (" << velocityModel.type() << ")" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

forAll(U, celli)
{
    U[celli] = velocityModel(mesh.C()[celli]);
}

forAll(U.boundaryField(), patchi)
{
    fvPatchVectorField& pU = U.boundaryFieldRef()[patchi];

    if (!pU.coupled())
    {
        forAll(pU, facei)
        {
            pU[facei] = velocityModel(pU.patch().Cf()[facei]);
        }
    }
}

U.correctBoundaryConditions();

// The face flux is evaluated directly from the analytical velocity at the face
// centres rather than interpolated from U. For the Taylor-Green vortex this
// makes the normal velocity on the domain boundary exactly zero, so the
// problem is closed without any boundary condition on the moments.
Info<< "Creating face flux phi" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimVolume/dimTime, 0)
);

forAll(phi, facei)
{
    phi[facei] = velocityModel(mesh.Cf()[facei]) & mesh.Sf()[facei];
}

forAll(phi.boundaryField(), patchi)
{
    fvsPatchScalarField& pPhi = phi.boundaryFieldRef()[patchi];
    const fvPatch& p = pPhi.patch();

    forAll(pPhi, facei)
    {
        pPhi[facei] = velocityModel(p.Cf()[facei]) & p.Sf()[facei];
    }
}

// The transported equation is conservative, so a non-solenoidal discrete flux
// is not an error - but for a field that is analytically divergence free it is
// a direct measure of the mesh resolution, and worth seeing.
{
    const volScalarField divPhi(fvc::div(phi));

    Info<< "    max |div(u)| = " << gMax(mag(divPhi.primitiveField()))
        << ", net flux through the domain = "
        << gSum(divPhi.primitiveField()*mesh.V()) << endl;
}

// * * * * * * * * * * * * * * *  Moment fields  * * * * * * * * * * * * * * //

// All moments share one set of boundary conditions and dimensions, so they are
// all built from a single template field 0/M rather than from nMoments
// near-identical files.
Info<< "\nCreating " << nMoments << " moment fields from the template 0/M"
    << endl;

volScalarField Mtemplate
(
    IOobject
    (
        "M",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

PtrList<volScalarField> M(nMoments);

forAll(M, k)
{
    M.set
    (
        k,
        new volScalarField
        (
            IOobject
            (
                "M" + Foam::name(k),
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            Mtemplate
        )
    );

    volScalarField& Mk = M[k];

    forAll(Mk, celli)
    {
        const point& p = mesh.C()[celli];

        scalar spatial = 1;

        if (spatialShape == "gaussianBlob")
        {
            spatial = std::exp
            (
                -magSqr(p - blobCentre)/(2*sigmaSpace*sigmaSpace)
            );
        }
        else if (spatialShape != "uniform")
        {
            FatalErrorInFunction
                << "Unknown spatialShape " << spatialShape << nl
                << "Valid shapes are: uniform, gaussianBlob"
                << exit(FatalError);
        }

        Mk[celli] = phaseMoments[k]*spatial;
    }

    Mk.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * *  Robustness  * * * * * * * * * * * * * * * //

const dictionary& repairDict = qmomProperties.subOrEmptyDict("repair");

// Last-resort repair of cells whose moment vector cannot be inverted at all -
// typically the "negative recurrence" failure, where m_0 is perfectly healthy
// but the sequence has lost realizability. Such a cell is rebuilt from the
// nominal phase-space shape rescaled to the zeroth moment it already carries,
// which is realizable by construction because the realizable set is a cone.
// This conserves m_0 exactly but overwrites the higher moments, which had
// stopped meaning anything anyway. Off by default: a run that needs it is
// telling you the resolution or the time step is marginal.
const Switch repairBrokenCells
(
    repairDict.lookupOrDefault<Switch>("repairBrokenCells", false)
);

const scalar m0Floor(repairDict.lookupOrDefault<scalar>("m0Floor", 0));

// Nominal shape: the initial distribution normalised to m_0 = 1
std::vector<double> nominalShape(phaseMoments);
for (label k = 0; k < nMoments; ++k)
{
    nominalShape[k] /= phaseMoments[0];
}

// * * * * * * * * * * * * * *  Analytical solution  * * * * * * * * * * * * //

const qmom::exactSolution exactSolution
(
    qmomProperties.subOrEmptyDict("errorAnalysis"),
    growthModelType,
    growthRate,
    phaseMoments,
    sigmaSpace,
    blobCentre,
    velocityModel
);

// Write the initial moment fields so that the initial condition is available
// to post-processing alongside the later times.
forAll(M, k)
{
    M[k].write();
}
