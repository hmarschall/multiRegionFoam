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

#include "fvCFD.H"
#include "interTrackFluid.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(interTrackFluid, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        interTrackFluid,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::interTrackFluid::interTrackFluid
(
    const Time& runTime,
    const word& regionName
)
:
    regionType(runTime, regionName),
    regionName_(regionName),

    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            mesh().time().constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    pimple_(mesh()),

    rho_(nullptr),
    mu_(nullptr),
    U_(nullptr),
    phi_(nullptr),
    p_(nullptr),

    pRefCell_(0),
    pRefValue_
    (
        readScalar(mesh().solutionDict().subDict("PISO").lookup("pRefValue"))
    ),

    mrfProperties_
    (
        IOobject
        (
            "mrfProperties",
            mesh().time().constant(),
            mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    movingReferenceFrame_
    (
        IOobject
        (
            "movingReferenceFrame",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    ),
    lambdaFf_
    (
        readScalar(mrfProperties_.lookup("lambdaFf"))
    ),
    lambdaF0_
    (
        readScalar(mrfProperties_.lookup("lambdaF0"))
    ),
    centerFromMesh_
    (
        mrfProperties_.lookupOrDefault<word>("centreFromMesh", regionName_)
    ),
    center_(vector::zero),
    center0_(vector::zero),
    XF_("XF", dimLength, vector::zero),
    UF_("UF", dimVelocity, vector::zero),
    aF_("aF", dimAcceleration, vector::zero),

    closedVolume_
    (
        mesh().solutionDict()
        .lookupOrDefault<Switch>("closedVolume", false)
    ),
    hasSpacePatch_
    (
        mesh().solutionDict()
        .lookupOrDefault<Switch>("hasSpacePatch", false)
    ),

    myTimeIndex_(mesh().time().timeIndex()),
    adjustTimeStep_
    (
        mesh().time().controlDict()
        .lookupOrDefault<Switch>("adjustTimeStep", false)
    ),
    maxCo_
    (
        mesh().time().controlDict().lookupOrDefault<scalar>("maxCo", 1.0)
    ),
    maxDeltaT_
    (
        mesh().time().controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT)
    ),
    tUEqn(),
    cumulativeContErr_(0)
{
    rho_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "rho",
        dimensionedScalar(transportProperties_.lookup("rho")),
        false
    );

    mu_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "mu",
        dimensionedScalar(transportProperties_.lookup("mu")),
        false
    );

    U_ = lookupOrRead<volVectorField>
    (
        mesh(),
        mesh().solutionDict().lookupOrDefault<word>("velocityName", "U"),
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

    p_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "p",
        true,
        true
    );

    IOdictionary fvSchemesDict
    (
        IOobject
        (
            "fvSchemes",
            mesh().time().system(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dictionary fluxRequiredDict = fvSchemesDict.subDict("fluxRequired");

    for(label i=0; i<fluxRequiredDict.size(); i++)
    {
        if(fluxRequiredDict.toc()[i] != "default")
        {
            mesh().schemesDict().setFluxRequired(fluxRequiredDict.toc()[i]);
        }
    }

    // Get pressure reference cell
#   include "setRefCell.H"

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::interTrackFluid::~interTrackFluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::interTrackFluid::correct()
{
    // do nothing, add as required
}


Foam::scalar Foam::regionTypes::interTrackFluid::getMinDeltaT()
{
    // do nothing, add as required
    return GREAT;
}


void Foam::regionTypes::interTrackFluid::setCoupledEqns()
{
    // do nothing, add as required
}

void Foam::regionTypes::interTrackFluid::postSolve()
{
#       include "updateMovingReferenceFrame.H"
}

void Foam::regionTypes::interTrackFluid::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::interTrackFluid::prePredictor()
{
    Info<< nl << "Pre-predictor for " << this->typeName
    << " in region " << mesh().name()
    << nl << endl;
}

void Foam::regionTypes::interTrackFluid::momentumPredictor()
{
    Info<< nl << "Momentum predictor for " << this->typeName
        << " in region " << mesh().name()
        << nl << endl;

    // Make the fluxes relative
    phi_() -= fvc::meshPhi(rho_(), U_());

    if (myTimeIndex_ < mesh().time().timeIndex())
    {
#       include "CourantNo.H"
    }

    tUEqn =
    (
        fvm::ddt(rho_(), U_())
      + rho_()*aF_
      + fvm::div(fvc::interpolate(rho_())*phi_(), U_(), "div(phi,U)")
      - fvm::laplacian(mu_(), U_())
    );

    fvVectorMatrix& UEqn = tUEqn();

    solve(UEqn == -fvc::grad(p_()));
}

void Foam::regionTypes::interTrackFluid::pressureCorrector()
{
    Info<< nl << "Pressure corrector for " << this->typeName
    << " in region " << mesh().name()
    << nl << endl;

    fvVectorMatrix& UEqn = tUEqn();

    // --- PISO loop
    while (pimple_.correct())
    {
        volScalarField AU = UEqn.A();

        U_() = UEqn.H()/AU;

        phi_() = (fvc::interpolate(U_()) & mesh().Sf());

#       include "correctPatchPhi.H"

        while (pimple_.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(1.0/AU, p_())
             == fvc::div(phi_())
            );

            // Find reference cell
            if (closedVolume_)
            {
                pEqn.setReference(pRefCell_, pRefValue_);
            }

            pEqn.solve();

            if (pimple_.finalNonOrthogonalIter())
            {
                phi_() -= pEqn.flux();
            }
        }

#               include "continuityErrs.H"

        // Momentum corrector
        U_() -= fvc::grad(p_())/AU;
        U_().correctBoundaryConditions();
    }

    myTimeIndex_ = mesh().time().timeIndex();

    Info<< nl
        << mesh().name() << " Pressure:" << nl
        << "  max: " << gMax(p_()) << nl
        << "  min: " << gMin(p_()) << nl
        << "  mean: " << gAverage(p_()) << nl
        << mesh().name() << " Velocity:" << nl
        << "  max: " << gMax(U_()) << nl
        << "  min: "<< gMin(U_()) << nl
        << "  mean: " << gAverage(U_()) << nl
        << mesh().name() << " Volume: "
        << gSum(mesh().V()) << nl
        << endl;
}

void Foam::regionTypes::interTrackFluid::meshMotionCorrector()
{
    mesh().update();

    label intPatchID = mesh().boundaryMesh().findPatchID("interface");

    if (intPatchID != -1)
    {
        Info<< "phi boundary field AFTER mesh motion : sum local = "
            << gSum(mag(phi_().boundaryField()[intPatchID]))
            << ", global = "
            << gSum(phi_().boundaryField()[intPatchID]) << endl;

        Info<< "mesh.phi boundary field AFTER mesh motion :"
            << " sum local = " << gSum(mag(mesh().phi().boundaryField()[intPatchID]))
            << ", global = " << gSum(mesh().phi().boundaryField()[intPatchID])
            << endl;

        Info<< "mesh.phi old time boundary field AFTER mesh motion :"
            << " sum local = " << gSum(mag(mesh().phi().oldTime().boundaryField()[intPatchID]))
            << ", global = " << gSum(mesh().phi().oldTime().boundaryField()[intPatchID])
            << endl;

        Info<< "fvc::meshPhi boundary field AFTER mesh motion : sum local = "
            << gSum(mag(fvc::meshPhi(rho_(),U_())().boundaryField()[intPatchID]))
            << ", global = "
            << gSum(fvc::meshPhi(rho_(),U_())().boundaryField()[intPatchID]) << endl;

        scalarField netPhiA =
            phi_().boundaryField()[intPatchID]
          - fvc::meshPhi(rho_(),U_())().boundaryField()[intPatchID];

        Info<< "Moving surface continuity error AFTER mesh motion : sum local = "
            << gSum(mag(netPhiA)) << ", global = " << gSum(netPhiA)
            << endl << endl;
    }
}

// ************************************************************************* //
