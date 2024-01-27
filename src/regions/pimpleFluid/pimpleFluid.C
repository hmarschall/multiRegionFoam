/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "label.H"
#include "pimpleFluid.H"

#include "fvCFD.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(pimpleFluid, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        pimpleFluid,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::pimpleFluid::pimpleFluid
(
    const Time& runTime,
    const word& regionName
)
:
    regionType(runTime, regionName),

    regionName_(regionName),

    pimple_(mesh()),

    laminarTransport_(nullptr),
    rho_(nullptr),

    turbulence_(nullptr),

    U_(nullptr),
    phi_(nullptr),
    pKin_(nullptr),
    p_(nullptr),

    sigma_(nullptr),

    rAU_(nullptr),

    pRefCell_
    (
        pimple_.dict().lookupOrDefault<label>("pRefCell", 0)
    ),
    pRefValue_
    (
        pimple_.dict().lookupOrDefault<scalar>("pRefValue", 0.0)
    ),

    mrfZones_(mesh()),
    myTimeIndex_(mesh().time().timeIndex()),

    maxCo_
    (
        mesh().time().controlDict().lookupOrDefault<scalar>("maxCo", 1.0)
    ),
    maxDeltaT_
    (
        mesh().time().controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT)
    ),

    correctPhi_
    (
        pimple_.dict().lookupOrDefault<Switch>("correctPhi", false)
    ),
    checkMeshCourantNo_
    (
        pimple_.dict().lookupOrDefault<Switch>("checkMeshCourantNo", false)
    ),

    sumLocalContErr_(0),
    globalContErr_(0),
    cumulativeContErr_(0)
{
    U_ = lookupOrRead<volVectorField>
    (
        mesh(),
        "U",
        true,
        true
    );

    pKin_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "pKin",
        true,
        true
    );

    phi_ = lookupOrRead<surfaceScalarField>
    (
        mesh(),
        "phi",
        false,
        true,
        linearInterpolate(U_()) & mesh().Sf()
    );

    laminarTransport_.set(new singlePhaseTransportModel(U_(), phi_()));

    rho_.set(new dimensionedScalar(laminarTransport_().lookup("rho")));

    turbulence_ = incompressible::turbulenceModel::New
    (
        U_(), phi_(), laminarTransport_()
    );

    p_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "p",
        false,
        true,
        rho_().value()*pKin_()
    );

    sigma_ = lookupOrRead<volSymmTensorField>
    (
        mesh(),
        "sigma",
        false,
        true,
        rho_().value()
       *(
            - pKin_()*symmTensor(1,0,0,1,0,1)
            - turbulence_().devReff()
        )
    );

    wordList rAUPatchFieldTypes
    (
        U_().boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );
    rAU_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "rAU",
        mesh().time().deltaT(),
        rAUPatchFieldTypes,
        true
    );

    setRefCell(pKin_(), pimple_.dict(), pRefCell_, pRefValue_);
    mesh().schemesDict().setFluxRequired(pKin_().name());
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::pimpleFluid::~pimpleFluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::pimpleFluid::correct()
{
#       include "pimpleFluidCourantNo.H"
}


Foam::scalar Foam::regionTypes::pimpleFluid::getMinDeltaT()
{
    //- TODO: implement deltaT based on CFL criteria
    return GREAT;
}


void Foam::regionTypes::pimpleFluid::setCoupledEqns()
{
    // do nothing, add as required
}


void Foam::regionTypes::pimpleFluid::postSolve()
{
    if (myTimeIndex_ < mesh().time().timeIndex())
    {
        mrfZones_.translationalMRFs().correctMRF();

        mrfZones_.translationalMRFs().correctBoundaryVelocity(U_(), phi_());

        myTimeIndex_ = mesh().time().timeIndex();
    }
}


void Foam::regionTypes::pimpleFluid::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::pimpleFluid::prePredictor()
{
    Info<< nl << "Pre-predictor for " << this->typeName
        << " in region " << mesh().name()
        << nl << endl;

    if (mesh().changing() && correctPhi_)
    {
#       include "pimpleFluidCorrectPhi.H"
    }
}

void Foam::regionTypes::pimpleFluid::momentumPredictor()
{
    Info<< nl << "Momentum predictor for " << this->typeName
        << " in region " << mesh().name()
        << nl << endl;

    // Time derivative matrix
    tddtUEqn = fvm::ddt(U_());
    fvVectorMatrix& ddtUEqn = tddtUEqn();

    // Convection-diffusion matrix
    tHUEqn =
        (
            fvm::div(phi_(), U_())
          + turbulence_().divDevReff()
        );
    fvVectorMatrix& HUEqn = tHUEqn();

    mrfZones_.translationalMRFs().addFrameAcceleration(ddtUEqn);

    if (pimple_.momentumPredictor())
    {
        solve(relax(ddtUEqn + HUEqn) == -fvc::grad(pKin_()));
    }
}

void Foam::regionTypes::pimpleFluid::pressureCorrector()
{
    Info<< nl << "Pressure corrector for " << this->typeName
        << " in region " << mesh().name()
        << nl << endl;

    // Get cached matricies from momentum predictor
    fvVectorMatrix& ddtUEqn = tddtUEqn();
    fvVectorMatrix& HUEqn = tHUEqn();

    // --- PISO loop
    while (pimple_.correct())
    {
        // Update pressure BCs
        pKin_().boundaryField().updateCoeffs();

        // Prepare clean 1/a_p without time derivative and under-relaxation
        // contribution
        rAU_() = 1.0/HUEqn.A();

        // Calculate U from convection-diffusion matrix
        U_() = rAU_()*HUEqn.H();

        // Consistently calculate flux
        pimple_.calcTransientConsistentFlux(phi_(), U_(), rAU_(), ddtUEqn);

        // Global flux balance
        adjustPhi(phi_(), U_(), pKin_());

        while (pimple_.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian
                (
                    fvc::interpolate(rAU_())/pimple_.aCoeff(U_().name()),
                    pKin_(),
                    "laplacian(rAU,pKin)"
                )
             ==
                fvc::div(phi_())
            );

            pEqn.setReference(pRefCell_, pRefValue_);
            pEqn.solve
            (
                mesh().solutionDict()
                .solver(pKin_().select(pimple_.finalInnerIter()))
            );

            if (pimple_.finalNonOrthogonalIter())
            {
                phi_() -= pEqn.flux();
            }
        }

        //- Pressure relaxation except for last corrector
        if (!pimple_.finalIter())
        {
            pKin_().relax();
        }

#       include "pimpleFluidMovingMeshContinuityErrs.H"

        // Consistently reconstruct velocity after pressure equation. Note: flux is
        // made relative inside the function
        pimple_.reconstructTransientVelocity(U_(), phi_(), ddtUEqn, rAU_(), pKin_());

        // Update pressure field
        p_() = rho_().value()*pKin_();

        // Update sigma field
        sigma_() = rho_().value()
           *(
                - pKin_()*symmTensor(1,0,0,1,0,1)
                - turbulence_().devReff()
            );
    }

    turbulence_().correct();

    Info<< nl
        << mesh().name() << " Pressure:" << nl
        << "  max: " << gMax(p_()) << nl
        << "  min: " << gMin(p_()) << nl
        << "  mean: " << gAverage(p_()) << nl
        << mesh().name() << " Velocity:" << nl
        << "  max: " << gMax(U_()) << nl
        << "  min: "<< gMax(U_()) << nl
        << "  mean: " << gAverage(U_()) << nl
        << mesh().name() << " Volume: "
        << gSum(mesh().V()) << nl
        << endl;
}

void Foam::regionTypes::pimpleFluid::meshMotionCorrector()
{
    // Make the fluxes absolute
    fvc::makeAbsolute(phi_(), U_());

    mesh().update();

#       include "pimpleFluidVolContinuity.H"

    if (mesh().changing() && correctPhi_)
    {
#       include "pimpleFluidCorrectPhi.H"
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi_(), U_());

    if (mesh().moving() && checkMeshCourantNo_)
    {
#           include "pimpleFluidMeshCourantNo.H"
    }

    if (mesh().changing())
    {
#           include "pimpleFluidCourantNo.H"
    }
}

// ************************************************************************* //
