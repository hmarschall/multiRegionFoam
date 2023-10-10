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

#include "icoFluid.H"

#include "fvCFD.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(icoFluid, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        icoFluid,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::icoFluid::icoFluid
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
    phiHbyA_(nullptr),
    p_(nullptr),
    UfHeader_
    (
        "Uf",
        mesh().time().timeName(),
        mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Uf_
    (
        UfHeader_,
        mesh(),
        dimensionedVector("0", dimVelocity, vector::zero)
    ),
    pcorrTypes_(),
    pcorr_(nullptr),
    sigma_(nullptr),
    UUrf_(1),

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
    pRefCell_(0),
    pRefValue_
    (
        readScalar(mesh().solutionDict().subDict("PISO").lookup("pRefValue"))
    ),
    innerResidual_(1),
    residualPressure_(1),

    corr_(0),
    corrNonOrtho_(0),

    mrfZones_(mesh()),
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
    )
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

    word velocityName
    (
        mesh().solutionDict()
        .lookupOrDefault<word>("velocityName", "U")
    );

    U_ = lookupOrRead<volVectorField>
    (
        mesh(),
        velocityName,
        true,
        true
    );

    UUrf_ = mesh().solutionDict().subDict("relaxationFactors")
            .lookupOrDefault<scalar>(velocityName, 1);

    phi_ = lookupOrRead<surfaceScalarField>
    (
        mesh(),
        "phi",
        false,
        true,
        linearInterpolate(U_()) & mesh().Sf()
    );
    phi_().oldTime();

    phiHbyA_ = lookupOrRead<surfaceScalarField>
    (
        mesh(),
        "phiHbyA",
        false,
        true,
        phi_()
    );

    p_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "p",
        true,
        true
    );

    pcorrTypes_ = wordList
    (
        p_().boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    for (label i = 0; i<p_().boundaryField().size(); i++)
    {
        if (p_().boundaryField()[i].fixesValue())
        {
            pcorrTypes_[i] = fixedValueFvPatchScalarField::typeName;
        }
    };

    pcorr_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "pcorr",
        dimensionedScalar("pcorr", p_().dimensions(), 0.0),
        pcorrTypes_,
        true
    );

    sigma_ = lookupOrRead<volSymmTensorField>
    (
        mesh(),
        "sigma",
        false,
        false,
        (-p_()*symmTensor(1,0,0,1,0,1))
      + (mu_()*twoSymm(fvc::grad(U_())))
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

Foam::regionTypes::icoFluid::~icoFluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::icoFluid::correct()
{
    // Get pressure reference cell
//#   include "setRefCell.H"

    if (mesh().changing())
    {
#       include "correctPhi.H"
    }
}


Foam::scalar Foam::regionTypes::icoFluid::getMinDeltaT()
{
    //- TODO: implement deltaT based on CFL criteria
    return GREAT;
}


void Foam::regionTypes::icoFluid::setCoupledEqns()
{
    // do nothing, add as required
}


void Foam::regionTypes::icoFluid::postSolve()
{
    if (myTimeIndex_ < mesh().time().timeIndex())
    {
        mrfZones_.translationalMRFs().correctMRF();

        mrfZones_.translationalMRFs().correctBoundaryVelocity(U_(), phi_());

        myTimeIndex_ = mesh().time().timeIndex();
    }
}


void Foam::regionTypes::icoFluid::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::icoFluid::prePredictor()
{
    Info<< nl << "Pre-predictor for " << this->typeName
        << " in region " << mesh().name()
        << nl << endl;

    if (myTimeIndex_ < mesh().time().timeIndex())
    {
#       include "CourantNo.H"
#       include "setDeltaT.H"
    }

    // Get pressure reference cell
//#   include "setRefCell.H"

    if (mesh().changing())
    {
#       include "correctPhi.H"
    }
}

void Foam::regionTypes::icoFluid::momentumPredictor()
{
    Info<< nl << "Momentum predictor for " << this->typeName
        << " in region " << mesh().name()
        << nl << endl;

    // Make the fluxes relative to the mesh motion
    phi_() == (phi_() - fvc::meshPhi(rho_(), U_()));

    // Convection-diffusion matrix
    fvVectorMatrix HUEqn
    (
        fvm::div(fvc::interpolate(rho_())*phi_(), U_(), "div(phi,U)")
      - fvm::laplacian(mu_(), U_())
    );

    // Time derivative matrix
    fvVectorMatrix ddtUEqn(fvm::ddt(rho_(), U_()));
    mrfZones_.translationalMRFs().addFrameAcceleration(ddtUEqn, rho_());

    // UEqn
    tUEqn = ddtUEqn + HUEqn;
    fvVectorMatrix& UEqn = tUEqn();

    // Save source and boundaryCoeffs
    vectorField S0 = UEqn.source();
    FieldField<Field, vector> B0 = UEqn.boundaryCoeffs();

    // Relax and solve momentum equation
    UEqn.relax(UUrf_);

    if (pimple_.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p_()));
    }

    // reset equation to ensure relaxation parameter
    // is not causing problems with time consistency
    UEqn = (ddtUEqn + HUEqn);
    UEqn.source() = S0;
    UEqn.boundaryCoeffs() = B0;
}

void Foam::regionTypes::icoFluid::pressureCorrector()
{
    Info<< nl << "Pressure corrector for " << this->typeName
        << " in region " << mesh().name()
        << nl << endl;

    fvVectorMatrix& UEqn = tUEqn();

    // --- PISO loop
    while (pimple_.correct())
    {
        p_().storePrevIter();

        volScalarField AU = UEqn.A();
        volVectorField HU = UEqn.H();

        U_() = HU/AU;

        phi_() =
        (
            (fvc::interpolate(HU)/fvc::interpolate(AU))
          & mesh().Sf()
        );

#       include "correctPatchPhi.H"

        while (pimple_.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(1.0/fvc::interpolate(AU), p_())
             == fvc::div(phi_())
            );

            // Find reference cell
            if (closedVolume_)
            {
                //#   include "setRefCell.H"
                pEqn.setReference(pRefCell_, pRefValue_);
            }

            pEqn.solve
            (
                 mesh().solutionDict()
                 .solver(p_().select(pimple_.finalInnerIter()))
            );

            if (pimple_.finalNonOrthogonalIter())
            {
                phi_() -= pEqn.flux();
            }
        }

        //- Pressure relaxation except for last corrector
        if (!pimple_.finalIter())
        {
            p_().relax();
        }

        // Momentum corrector
        U_() -= fvc::grad(p_())/AU;

        U_().correctBoundaryConditions();
        p_().correctBoundaryConditions();

        // Update sigma field
        sigma_() =
            (-p_()*symmTensor(1,0,0,1,0,1))
          + (mu_()*twoSymm(fvc::grad(U_())));
    }

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

void Foam::regionTypes::icoFluid::meshMotionCorrector()
{
    mesh().update();

    if(!UfHeader_.headerOk())
    {
        Uf_ = fvc::interpolate(U_());
        surfaceVectorField n(mesh().Sf()/mesh().magSf());
        Uf_ += n*(phi_()/mesh().magSf() - (n & Uf_));
    }

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
