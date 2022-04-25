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
#include "navierStokesFluid.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "pimpleControl.H"
#include "solutionControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(navierStokesFluid, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        navierStokesFluid,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::navierStokesFluid::navierStokesFluid
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
            runTime.constant(),
            runTime, 
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    rhoFluid_
    (
        transportProperties_.subDict(regionName_).lookup("rho")
    ),
    rho_
    (
        IOobject
        (
            "rho",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        rhoFluid_
    ),
    muFluid_
    (
        transportProperties_.subDict(regionName_).lookup("mu")
    ),
    mu_
    (
        IOobject
        (
            "mu",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        muFluid_
    ),
    U_
    (
        IOobject
        (
           "U",
            mesh().time().timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    phi_
    (
		IOobject
		(
			"phi",
            mesh().time().timeName(),
            mesh(),
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		linearInterpolate(U_) & mesh().Sf()    
    ),
    phiHbyA_
    (
		IOobject
		(
			"phiHbyA",
            mesh().time().timeName(),
            mesh(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		phi_   
    ),
    p_
    (
		IOobject
		(
			"p",
            mesh().time().timeName(),
            mesh(),
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh()
    ),    
    AU_
    (
        IOobject
        (
            "AU",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimMass/(dimVolume*dimTime),
        zeroGradientFvPatchScalarField::typeName
    ),
    HU_
    (
        IOobject
        (
            "HU",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimMass/sqr(dimLength*dimTime),
        zeroGradientFvPatchVectorField::typeName
    ),              
    gradp_
    (
		IOobject
		(
			"grad(p)",
            mesh().time().timeName(),
            mesh(),
			IOobject::NO_READ,
            IOobject::NO_WRITE
		),
		fvc::grad(p_)   
    ), 
    gradU_
    (
		IOobject
		(
			"grad(U)",
            mesh().time().timeName(),
            mesh(),
			IOobject::NO_READ,
            IOobject::NO_WRITE
		),
		fvc::grad(U_) 
    ), 
    pcorrTypes_
    (
        p_.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    ),   
    pcorr_
    (
        IOobject
        (
            "pcorr",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("pcorr", p_.dimensions(), 0.0),
        pcorrTypes_
    ),
    UUrf_ 
    ( 
        mesh().solutionDict().subDict("relaxationFactors")
        .lookupOrDefault<scalar>(U_.name(), 1)   
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
    
    for (label i = 0; i<p_.boundaryField().size(); i++)
    {
        if (p_.boundaryField()[i].fixesValue())
        {
            pcorrTypes_[i] = fixedValueFvPatchScalarField::typeName;
        }
    };  
    
    gradU_.checkIn();
    gradp_.checkIn();

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
}

// left from createFields    
//#   include "createUf.H"
//#   include "createSf.H"
//#   include "setRefCell.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::navierStokesFluid::~navierStokesFluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::navierStokesFluid::correct()
{
    // do nothing, add as required
}


void Foam::regionTypes::navierStokesFluid::setRDeltaT()
{
    // do nothing, add as required
}


void Foam::regionTypes::navierStokesFluid::setCoupledEqns()
{  
    // do nothing, add as required
}


void Foam::regionTypes::navierStokesFluid::updateFields()
{
    // do nothing, add as required
}


void Foam::regionTypes::navierStokesFluid::solveRegion()
{
    Info << nl << "Solving for " << mesh().name() << endl;

    // --- PIMPLE loop      
    pimpleControl pimple(mesh());

    if (myTimeIndex_ < mesh().time().timeIndex())
    {
#       include "CourantNo.H"
#       include "setDeltaT.H"
    }

    // Get pressure reference cell
#   include "setRefCell.H"

    if (mesh().changing())
    {
#       include "correctPhi.H"
    }

    if (myTimeIndex_ < mesh().time().timeIndex())
    {
        mrfZones_.translationalMRFs().correctMRF();
    }

    while (pimple.loop())
    {
        // Make the fluxes relative to the mesh motion
        phi_ == (phi_ - fvc::meshPhi(rho_, U_));
//        fvc::makeRelative(phi_, U_);

//        mrfZones_.translationalMRFs().correctBoundaryVelocity(U_, phi_);

        U_.storePrevIter();

        // Convection-diffusion matrix
        fvVectorMatrix HUEqn
        (
            fvm::div(fvc::interpolate(rho_)*phi_, U_, "div(phi,U)")
          - fvm::laplacian(mu_, U_)
        );

        // Time derivative matrix
        fvVectorMatrix ddtUEqn(fvm::ddt(rho_, U_));

//        if (myTimeIndex_ < mesh().time().timeIndex())
        {
            mrfZones_.translationalMRFs().addFrameAcceleration(ddtUEqn, rho_);

//            myTimeIndex_ = mesh().time().timeIndex();
        }

        // UEqn
        fvVectorMatrix UEqn(ddtUEqn + HUEqn);

        // Save source and boundaryCoeffs
        vectorField S0 = UEqn.source();
        FieldField<Field, vector> B0 = UEqn.boundaryCoeffs();


        // Relax and solve momentum equation
        UEqn.relax(UUrf_);

        if (pimple.momentumPredictor()) 
        {
            solve(UEqn == -gradp_);
        }
        
        // reset equation to ensure relaxation parameter 
        // is not causing problems with time consistency
        UEqn = (ddtUEqn + HUEqn);
        UEqn.source() = S0;
        UEqn.boundaryCoeffs() = B0;


        // --- PISO loop
        while (pimple.correct())
        {
            p_.storePrevIter();   
                 
            AU_ = UEqn.A();
            HU_ = UEqn.H();

            U_ = HU_/AU_;

            phi_ =
            (
                (fvc::interpolate(HU_)/fvc::interpolate(AU_))
               & mesh().Sf()
            );

//            phi_ += fvc::ddtPhiCorr((1.0/AU_)(), rho_, U_, phi_);

//            phi += fvc::ddtPhiCorr
//                (
//                    (rho_/AU_)(),
//                    U_,
//                    phi_
//                );

            // Global flux continuity
//            adjustPhi(phiHbyA_, U_, p_);

            if (closedVolume_)
            {
                label intPatchID_ = 
                    mesh().boundaryMesh().findPatchID("interfaceShadow");

                forAll(phi_.boundaryField(), patchI)
                {
                    if
                    (
                        !phi_.boundaryField()[patchI].coupled()
                     && patchI != intPatchID_
//                     && isA<zeroGradientFvPatchVectorField>
//                        (
//                            U_.boundaryField()[patchI]
//                        )
                    )
                    {
                        phi_.boundaryField()[patchI] ==
                        (
                            U_.boundaryField()[patchI]
                            & mesh().Sf().boundaryField()[patchI]
                        );
                    }
                }
            }

            if (!closedVolume_)
            {
                forAll(phi_.boundaryField(), patchI)
                {
                    if
                    (
                        !phi_.boundaryField()[patchI].coupled()
                     && !p_.boundaryField()[patchI].fixesValue()
                    )
                    {
                        phi_.boundaryField()[patchI] ==
                        (
                            U_.boundaryField()[patchI]
                            & mesh().Sf().boundaryField()[patchI]
                        );
                    }
                }
            }

            if (closedVolume_)
            {
                label intPatchID_ = 
                    mesh().boundaryMesh().findPatchID("interfaceShadow");

                phi_.boundaryField()[intPatchID_] =
                    (
                        U_.boundaryField()[intPatchID_]
                      & mesh().Sf().boundaryField()[intPatchID_]
                    );

                scalarField weights =
                    mag(phi_.boundaryField()[intPatchID_]);

                if(mag(gSum(weights)) > VSMALL)
                {
                    weights /= gSum(weights);
                }

                phi_.boundaryField()[intPatchID_] -=
                    weights*gSum(phi_.boundaryField()[intPatchID_]);

                phi_.boundaryField()[intPatchID_] +=
                    p_.boundaryField()[intPatchID_].snGrad()
                   *mesh().magSf().boundaryField()[intPatchID_]
                   /AU_.boundaryField()[intPatchID_];
            }

            if (!closedVolume_ && hasSpacePatch_)
            {
                // get space patch index
                label scalePatchID =
                    mesh().boundaryMesh().findPatchID("space");

                //- Non-const access to flux on patch
                fvsPatchField<scalar>& phip = 
                    const_cast<fvsPatchField<scalar>& >
                    (
                        mesh().lookupObject<surfaceScalarField>("phi")
                        .boundaryField()[scalePatchID]
                    );

                scalar inletFlux = gSum(neg(phip)*phip);

                scalar outletFlux = gSum(pos(phip)*phip);

                if(outletFlux < VSMALL)
                {
                    outletFlux = VSMALL;
                }

                scalar outflowScaling = -inletFlux/outletFlux;

                phip += pos(phip)*phip*(outflowScaling - 1.0);
            }

            while (pimple.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(1.0/fvc::interpolate(AU_), p_) 
                 == fvc::div(phi_)
                );

               if (closedVolume_)
               {
                   pEqn.setReference(pRefCell_, pRefValue_);
               }

                pEqn.solve
                (
                    mesh().solutionDict()
                    .solver(p_.select(pimple.finalInnerIter()))
                );

                if (pimple.finalNonOrthogonalIter())
                {
                    phi_ -= pEqn.flux();
                }                            
            }

            //- Pressure relaxation except for last corrector
            if (!pimple.finalIter())
            {
                p_.relax();
            }

            // Update of pressure gradient
            gradp_ = fvc::grad(p_);

            // Momentum corrector
            U_ -= fvc::grad(p_)/AU_;

            U_.correctBoundaryConditions();

            // Update of velocity gradient
            gradU_ = fvc::grad(U_);       
        }

        {
            mrfZones_.translationalMRFs().correctBoundaryVelocity(U_, phi_);

            myTimeIndex_ = mesh().time().timeIndex();
        }
    }
}

// ************************************************************************* //
