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
    nu_(nullptr),
    U_(nullptr),
    phi_(nullptr),
    p_(nullptr),

    rAU_
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        runTime.deltaT(),
        zeroGradientFvPatchScalarField::typeName
    ),

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

    nu_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "nu", 
        dimensionedScalar(transportProperties_.lookup("nu")),
        false
    );

    U_ = lookupOrRead<volVectorField>
    (
        mesh(), 
        "U"
    );

    phi_ = lookupOrRead<surfaceScalarField>
    (
        mesh(),
        "phi",
        false,
        true,
        linearInterpolate(U_()) & mesh().Sf()
    );

    phi_().oldTime();

    // look up pressure field from object registry
    p_ = lookupOrRead<volScalarField>
    (
        mesh(), 
        "p"
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

    // left from createFields    
//#   include "createUf.H"
//#   include "createSf.H"
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
    // do nothing, add as required
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

    if (myTimeIndex_ < mesh().time().timeIndex())
    {
        mrfZones_.translationalMRFs().correctMRF();
    }
}

void Foam::regionTypes::icoFluid::momentumPredictor()
{
    Info<< nl << "Momentum predictor for " << this->typeName
        << " in region " << mesh().name() 
        << nl << endl;

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi_(), U_());

        U_().storePrevIter();

        // Time derivative matrix
        tddtUEqn_ = fvm::ddt(U_());
        fvVectorMatrix& ddtUEqn = tddtUEqn_();

        mrfZones_.translationalMRFs().addFrameAcceleration(ddtUEqn);

        // Convection-diffusion matrix
        tHUEqn_ =
        (
            fvm::div(phi_(), U_(), "div(phi,U)")
          - fvm::laplacian(nu_(), U_())
        );
        fvVectorMatrix& HUEqn = tHUEqn_();

        if (pimple_.momentumPredictor()) 
        {
            solve(ddtUEqn + HUEqn == -fvc::grad(p_()));
        }

        // Store rAU for pressure correction 
        //rAU_ = 1.0/HUEqn.A();

        // Prepare clean 1/a_p without time derivative contribution
        rAU_ = 1.0/HUEqn.A();
}

void Foam::regionTypes::icoFluid::pressureCorrector()
{
    Info<< nl << "Pressure corrector for " << this->typeName
        << " in region " << mesh().name() 
        << nl << endl;

    fvVectorMatrix& ddtUEqn = tddtUEqn_();
    fvVectorMatrix& HUEqn = tHUEqn_();

    // --- PISO loop
    while (pimple_.correct())
    {
        p_().storePrevIter();   

        // Calculate U from convection-diffusion matrix
        U_() = rAU_*HUEqn.H();

        pimple_.calcTransientConsistentFlux(phi_(), U_(), rAU_, ddtUEqn);
        // phi_() =
        // (
        //     fvc::interpolate(U_()) & mesh().Sf()
        // );

#       include "correctPatchPhi.H"

        while (pimple_.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian
                (
                    fvc::interpolate(rAU_)/pimple_.aCoeff("U"),
                    p_(),
                    "laplacian(rAU,p)"
                )
             == fvc::div(phi_())
            );

            // fvScalarMatrix pEqn
            // (
            //     fvm::laplacian(1.0/fvc::interpolate(AU), p_()) 
            //  == fvc::div(phi_())
            // );

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
        //U_() -= fvc::grad(p_())*rAU_;
        pimple_.reconstructTransientVelocity(U_(), phi_(), ddtUEqn, rAU_, p_());

        U_().correctBoundaryConditions();
        p_().correctBoundaryConditions();  
    }

    // set mrf boundary velocity
    mrfZones_.translationalMRFs().correctBoundaryVelocity(U_(), phi_());

    // set local time index after mrf boundary velocity has been set 
    myTimeIndex_ = mesh().time().timeIndex();

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

    // if(!UfHeader_.headerOk())
    // {
    //     Uf_ = fvc::interpolate(U_());
    //     surfaceVectorField n(mesh().Sf()/mesh().magSf());
    //     Uf_ += n*(phi_()/mesh().magSf() - (n & Uf_));
    // }
}

// ************************************************************************* //
