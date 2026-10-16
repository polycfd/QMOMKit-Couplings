/*---------------------------------------------------------------------------*\
Class
    Foam::qmom models

Description
    Growth laws, velocity fields and analytical solutions for
    qmomAdvectionFoam. Each is selected by name from constant/qmomProperties:
    growthModel (linear, constant, none), velocityModel (noFlow, constantX,
    taylorGreen) and errorAnalysis (none, constantSpace, advectedGaussian,
    swirlGaussian).

\*---------------------------------------------------------------------------*/

#ifndef qmomModels_H
#define qmomModels_H

#include "dictionary.H"
#include "vector.H"
#include "mathematicalConstants.H"
#include "error.H"

#include <qmomkit/qmomkit.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

namespace Foam
{
namespace qmom
{

// * * * * * * * * * * * * * *  Phase-space growth  * * * * * * * * * * * * * //

//- Build the phase-space velocity G(t, xi) = dxi/dt handed to QMOMKit.
//
//  The population balance equation solved in phase space is
//
//      df/dt + d/dxi( G(xi) f ) = 0
//
//  from which QMOMKit builds the moment rates
//
//      dM_k/dt = sum_i w_i k G(xi_i) xi_i^(k-1)
//
//  once it has inverted the moments into a quadrature.
//
//  QMOMKit evaluates G for all quadrature nodes in one call: xi holds the
//  node coordinates packed node-major and G, already sized by the library,
//  receives dxi/dt in the same layout. Every law below acts on each
//  coordinate independently, so the packing does not matter here.
//
//  linear   G(xi) = growthRate*xi -> dM_k/dt = k*growthRate*M_k exactly, so
//                   M_k(t) = M_k(0) exp(k*growthRate*t). The closure is exact
//                   for this law, which is what makes problem A a clean test
//                   of the QMOM machinery on its own.
//  constant G(xi) = growthRate    -> dM_k/dt = k*growthRate*M_(k-1). The
//                   distribution is shifted rigidly, f(xi, t) = f(xi - g t, 0),
//                   and the closure is exact for this law too.
//  none     G(xi) = 0             -> the moments are untouched in phase space.
inline qmomkit::PhaseSpaceVelocity newGrowthModel
(
    const word& type,
    const scalar growthRate
)
{
    const double g = growthRate;

    if (type == "linear")
    {
        return [g]
        (
            double,
            const std::vector<double>& xi,
            std::vector<double>& G,
            std::size_t
        )
        {
            for (std::size_t i = 0; i < xi.size(); ++i)
            {
                G[i] = g*xi[i];
            }
        };
    }
    else if (type == "constant")
    {
        return [g]
        (
            double,
            const std::vector<double>&,
            std::vector<double>& G,
            std::size_t
        )
        {
            std::fill(G.begin(), G.end(), g);
        };
    }
    else if (type == "none")
    {
        return []
        (
            double,
            const std::vector<double>&,
            std::vector<double>& G,
            std::size_t
        )
        {
            std::fill(G.begin(), G.end(), 0.0);
        };
    }

    FatalErrorInFunction
        << "Unknown growthModel " << type << nl
        << "Valid types are: linear, constant, none"
        << exit(FatalError);

    return nullptr;
}


// * * * * * * * * * * * * * *  Velocity models  * * * * * * * * * * * * * * //

//- Steady physical-space velocity field u(x).
class velocityModel
{
    word type_;
    scalar U0_;
    scalar omega_;

public:

    explicit velocityModel(const dictionary& dict)
    :
        type_(dict.lookup<word>("type")),
        U0_(dict.lookupOrDefault<scalar>("U0", 1)),
        omega_(dict.lookupOrDefault<scalar>("omega", 1))
    {
        if (type_ != "noFlow" && type_ != "constantX" && type_ != "taylorGreen")
        {
            FatalIOErrorInFunction(dict)
                << "Unknown velocity model " << type_ << nl
                << "Valid types are: noFlow, constantX, taylorGreen"
                << exit(FatalIOError);
        }
    }

    const word& type() const
    {
        return type_;
    }

    scalar U0() const
    {
        return U0_;
    }

    //- u(x). The Taylor-Green vortex uses k = pi/2 on [-1,1]^2, which puts a
    //  velocity node on every domain boundary, so the normal velocity there
    //  is exactly zero and no boundary treatment is needed.
    vector operator()(const point& p) const
    {
        if (type_ == "constantX")
        {
            return vector(U0_, 0, 0);
        }
        else if (type_ == "taylorGreen")
        {
            const scalar k = 0.5*constant::mathematical::pi;
            return vector
            (
               -omega_*std::cos(k*p.x())*std::sin(k*p.y()),
                omega_*std::sin(k*p.x())*std::cos(k*p.y()),
                0
            );
        }

        return Zero;
    }
};


// * * * * * * * * * * * * *  Analytical solutions  * * * * * * * * * * * * * //

//- Analytical raw moments M_k(x, t) of the four reference problems.
//
//  Every problem starts from X ~ N(mu, sigma) in phase space (QMOMKit
//  supplies E[X^k]) and the growth law does not depend on x, so M_k is the
//  product of a phase-space factor, given by phaseMoment(), and a spatial
//  factor. The types differ only in the spatial factor:
//
//      uniform          1                     (problem A)
//      advectedGaussian a Gaussian blob rigidly translated by u = U0 e_x and
//                       wrapped into the periodic box  (problems B and C)
//      swirlGaussian    the initial blob read at the foot of the characteristic
//                       of a steady velocity field     (problem D)
class exactSolution
{
    word type_;
    word growthType_;
    scalar growthRate_;

    //- E[X^k] of the initial phase-space distribution, from QMOMKit
    std::vector<double> phaseMoments_;

    scalar sigmaSpace_;
    point centre_;
    scalar length_;
    scalar traceStepSize_;
    label reportEvery_;

    const velocityModel& U_;

public:

    exactSolution
    (
        const dictionary& dict,
        const word& growthType,
        const scalar growthRate,
        const std::vector<double>& phaseMoments,
        const scalar sigmaSpace,
        const point& centre,
        const velocityModel& U
    )
    :
        type_(dict.lookupOrDefault<word>("type", "none")),
        growthType_(growthType),
        growthRate_(growthRate),
        phaseMoments_(phaseMoments),
        sigmaSpace_(sigmaSpace),
        centre_(centre),
        length_(dict.lookupOrDefault<scalar>("periodicLength", 2)),
        traceStepSize_(dict.lookupOrDefault<scalar>("traceStepSize", 0.01)),
        reportEvery_(dict.lookupOrDefault<label>("reportEvery", 1)),
        U_(U)
    {
        if
        (
            type_ != "none" && type_ != "constantSpace"
         && type_ != "advectedGaussian" && type_ != "swirlGaussian"
        )
        {
            FatalIOErrorInFunction(dict)
                << "Unknown exact solution " << type_ << nl
                << "Valid types are: none, constantSpace, advectedGaussian, "
                << "swirlGaussian" << exit(FatalIOError);
        }
    }

    bool active() const
    {
        return type_ != "none";
    }

    const word& type() const
    {
        return type_;
    }

    //- How often, in time steps, the error is worth reporting. The swirl
    //  reference costs a backward trace per cell and its cost grows linearly
    //  with t, so problem D reports far less often than every step.
    label reportEvery() const
    {
        return max(reportEvery_, label(1));
    }

    //- Trace a point back along the characteristic of the steady velocity
    //  field to where it was at t = 0, by integrating dX/ds = -u(X) with
    //  RK4. Taylor-Green trajectories have no elementary closed form (they
    //  reduce to elliptic integrals), so problem D's reference has to be
    //  computed this way. Valid only because the field is steady.
    point traceBackToOrigin(const point& p, const scalar t) const
    {
        if (t <= 0)
        {
            return p;
        }

        // The number of steps grows with t so that the step size, and hence
        // the accuracy of the trace, is the same at every reported time.
        const label nSteps = max(label(std::ceil(t/traceStepSize_)), label(1));
        const scalar h = t/nSteps;
        point x = p;

        for (label i = 0; i < nSteps; ++i)
        {
            const vector k1 = -U_(x);
            const vector k2 = -U_(x + 0.5*h*k1);
            const vector k3 = -U_(x + 0.5*h*k2);
            const vector k4 = -U_(x + h*k3);

            x += (h/6.0)*(k1 + 2*k2 + 2*k3 + k4);
        }

        return x;
    }

    //- Spatial factor of moment k at (x, t)
    scalar spatialFactor(const point& p, const scalar t) const
    {
        if (type_ == "constantSpace")
        {
            return 1;
        }
        else if (type_ == "advectedGaussian")
        {
            // Rigid translation by u = U0 e_x, wrapped into the periodic box
            point q = p - centre_;
            q.x() -= U_.U0()*t;
            q.x() -= length_*std::floor((q.x() + 0.5*length_)/length_);

            return std::exp(-magSqr(q)/(2*sigmaSpace_*sigmaSpace_));
        }
        else if (type_ == "swirlGaussian")
        {
            // The blob is the initial field read at the foot of the
            // characteristic. Since div(u) = 0 for the Taylor-Green vortex the
            // transport reduces to dM_k/dt = k*rate*M_k along a trajectory, so
            // composing the initial field with the flow map is exact. The blob
            // shears and does not stay Gaussian.
            const point q = traceBackToOrigin(p, t) - centre_;

            return std::exp(-magSqr(q)/(2*sigmaSpace_*sigmaSpace_));
        }

        return 1;
    }

    //- Phase-space factor of moment k at time t, i.e. E[X(t)^k]
    //
    //  linear   X(t) = X e^(rate t), so E[X^k] e^(k rate t).
    //  constant X(t) = X + rate t, a rigid shift, so by the binomial theorem
    //           sum_j C(k,j) E[X^j] (rate t)^(k-j). Valid for any initial
    //           distribution, not only the Gaussian.
    //  none     E[X^k], unchanged.
    scalar phaseMoment(const label k, const scalar t) const
    {
        if (growthType_ == "linear")
        {
            return phaseMoments_[k]*std::exp(scalar(k)*growthRate_*t);
        }
        else if (growthType_ == "constant")
        {
            const scalar shift = growthRate_*t;

            // Walk j from k down to 0 so that the binomial coefficient
            // C(k,j) and the power shift^(k-j) are both built incrementally
            scalar sum = 0;
            scalar binomial = 1;
            scalar shiftPower = 1;

            for (label j = k; j >= 0; --j)
            {
                sum += binomial*phaseMoments_[j]*shiftPower;

                binomial *= scalar(j)/scalar(k - j + 1);
                shiftPower *= shift;
            }

            return sum;
        }

        return phaseMoments_[k];
    }

    //- M_k(x, t)
    scalar value(const label k, const point& p, const scalar t) const
    {
        return phaseMoment(k, t)*spatialFactor(p, t);
    }
};

} // End namespace qmom
} // End namespace Foam

#endif

// ************************************************************************* //
