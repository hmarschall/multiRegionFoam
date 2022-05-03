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
    rho_(nullptr),
    muFluid_
    (
        transportProperties_.subDict(regionName_).lookup("mu")
    ),
    mu_(nullptr),
    velocityName_
    (
        mesh().solutionDict()
        .lookupOrDefault<word>("velocityName", "U")
    ),
    U_(nullptr),
    phi_(nullptr),
    phiHbyA_(nullptr),
    p_(nullptr),    
    AU_(nullptr),
    HU_(nullptr),              
    gradp_(nullptr), 
    gradU_(nullptr), 
    pcorrTypes_(),   
    pcorr_(nullptr),
    UUrf_ 
    ( 
        mesh().solutionDict().subDict("relaxationFactors")
        .lookupOrDefault<scalar>(velocityName_, 1)   
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
    pRefValue_(0),     
    
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
    // look up desity field from object registry
    if (mesh().foundObject<volScalarField>("rho"))
    {
        Info << nl << "Using already existing desity field in region "
             << mesh().name()
             << " for regionType "
             << this->name()
             << nl << endl;

        rho_.reset
        (
            const_cast<volScalarField*>
            (
                &mesh().lookupObject<volScalarField>("rho")
            )
        );
    }
    else // create new desity field
    {
        rho_.reset
        (
            new volScalarField
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
            )
        );
    }

    // look up viscosity field from object registry
    if (mesh().foundObject<volScalarField>("mu"))
    {
        Info << nl << "Using already existing viscosity field in region "
             << mesh().name()
             << " for regionType "
             << this->name()
             << nl << endl;

        mu_.reset
        (
            const_cast<volScalarField*>
            (
                &mesh().lookupObject<volScalarField>("mu")
            )
        );
    }
    else // create new desity field
    {
        mu_.reset
        (
            new volScalarField
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
            )
        );
    }

    // look up velocity field from object registry
    if (mesh().foundObject<volVectorField>("U"))
    {
        Info << nl << "Using already existing velocity field in region "
             << mesh().name()
             << " for regionType "
             << this->name()
             << nl << endl;

        U_.reset
        (
            const_cast<volVectorField*>
            (
                &mesh().lookupObject<volVectorField>("U")
            )
        );
    }
    else // read velocity field
    {
        U_.reset
        (
            new volVectorField
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
            )
        );
    }

    // look up flux field from object registry
    if (mesh().foundObject<surfaceScalarField>("phi"))
    {
        Info << nl << "Using already existing flux field in region "
             << mesh().name()
             << " for regionType "
             << this->name()
             << nl << endl;

        phi_.reset
        (
            const_cast<surfaceScalarField*>
            (
                &mesh().lookupObject<surfaceScalarField>("phi")
            )
        );
    }
    else // use pre-set velocity field
    {
        phi_.reset
        (
            new surfaceScalarField
            (
		        IOobject
		        (
			        "phi",
                    mesh().time().timeName(),
                    mesh(),
			        IOobject::NO_READ,
			        IOobject::AUTO_WRITE
		        ),
		        linearInterpolate(U_()) & mesh().Sf()
            )
        );
    }

    // look up phiHbyA_ field from object registry
    if (mesh().foundObject<surfaceScalarField>("phiHbyA"))
    {
        Info << nl << "Using already existing phiHbyA field in region "
             << mesh().name()
             << " for regionType "
             << this->name()
             << nl << endl;

        phiHbyA_.reset
        (
            const_cast<surfaceScalarField*>
            (
                &mesh().lookupObject<surfaceScalarField>("phiHbyA")
            )
        );
    }
    else // use pre-set flux field
    {
        phiHbyA_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "phiHbyA",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phi_()
            )
        );
    }

    // look up pressure field from object registry
    if (mesh().foundObject<volScalarField>("p"))
    {
        Info << nl << "Using already existing pressure field in region "
             << mesh().name()
             << " for regionType "
             << this->name()
             << nl << endl;

        p_.reset
        (
            const_cast<volScalarField*>
            (
                &mesh().lookupObject<volScalarField>("p")
            )
        );
    }
    else // create new pressure field
    {
        p_.reset
        (
            new volScalarField
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
            )
        );
    }

    // look up AU field from object registry
    if (mesh().foundObject<volScalarField>("AU"))
    {
        Info << nl << "Using already existing AU field in region "
             << mesh().name()
             << " for regionType "
             << this->name()
             << nl << endl;

        AU_.reset
        (
            const_cast<volScalarField*>
            (
                &mesh().lookupObject<volScalarField>("AU")
            )
        );
    }
    else // create new AU field
    {
        AU_.reset
        (
            new volScalarField
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
            )
        );
    }
 
    // look up HU field from object registry
    if (mesh().foundObject<volVectorField>("HU"))
    {
        Info << nl << "Using already existing HU field in region "
             << mesh().name()
             << " for regionType "
             << this->name()
             << nl << endl;

        HU_.reset
        (
            const_cast<volVectorField*>
            (
                &mesh().lookupObject<volVectorField>("HU")
            )
        );
    }
    else // create new HU field
    {
        HU_.reset
        (
            new volVectorField
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
            )
        );
    }

    // look up pressure gradient field from object registry
    if (mesh().foundObject<volVectorField>("grad(p)"))
    {
        Info << nl << "Using already existing pressure gradient field in region "
             << mesh().name()
             << " for regionType "
             << this->name()
             << nl << endl;

        gradp_.reset
        (
            const_cast<volVectorField*>
            (
                &mesh().lookupObject<volVectorField>("grad(p)")
            )
        );
    }
    else // create new pressure gradient field
    {
        gradp_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "grad(p)",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::grad(p_())   
            )
        );
    }

    // look up velocity gradient field from object registry
    if (mesh().foundObject<volTensorField>("grad(U)"))
    {
        Info << nl << "Using already existing velocity gradient field in region "
             << mesh().name()
             << " for regionType "
             << this->name()
             << nl << endl;

        gradU_.reset
        (
            const_cast<volTensorField*>
            (
                &mesh().lookupObject<volTensorField>("grad(U)")
            )
        );
    }
    else // create new velocity gradient field
    {
        gradU_.reset
        (
            new volTensorField
            (
                IOobject
                (
                    "grad(U)",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::grad(U_()) 
            )
        );
    }

    pcorrTypes_ = wordList
                (
                    p_().boundaryField().size(),
                    zeroGradientFvPatchScalarField::typeName
                );

    // look up pressure correction field from object registry
    if (mesh().foundObject<volScalarField>("pcorr"))
    {
        Info << nl << "Using already existing pressure correction field in region "
             << mesh().name()
             << " for regionType "
             << this->name()
             << nl << endl;

        pcorr_.reset
        (
            const_cast<volScalarField*>
            (
                &mesh().lookupObject<volScalarField>("pcorr")
            )
        );
    }
    else // create new pressure correction field
    {
        pcorr_.reset
        (
            new volScalarField
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
                dimensionedScalar("pcorr", p_().dimensions(), 0.0),
                pcorrTypes_
            )
        );
    }

    for (label i = 0; i<p_().boundaryField().size(); i++)
    {
        if (p_().boundaryField()[i].fixesValue())
        {
            pcorrTypes_[i] = fixedValueFvPatchScalarField::typeName;
        }
    };  
    
    gradU_().checkIn();
    gradp_().checkIn();

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
        phi_() == (phi_() - fvc::meshPhi(rho_(), U_()));
//        fvc::makeRelative(phi_(), U_());

//        mrfZones_.translationalMRFs().correctBoundaryVelocity(U_(), phi_());

        U_().storePrevIter();

        // Convection-diffusion matrix
        fvVectorMatrix HUEqn
        (
            fvm::div(fvc::interpolate(rho_())*phi_(), U_(), "div(phi,U)")
          - fvm::laplacian(mu_(), U_())
        );

        // Time derivative matrix
        fvVectorMatrix ddtUEqn(fvm::ddt(rho_(), U_()));

//        if (myTimeIndex_ < mesh().time().timeIndex())
        {
            mrfZones_.translationalMRFs().addFrameAcceleration(ddtUEqn, rho_());

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
            solve(UEqn == -gradp_());
        }
        
        // reset equation to ensure relaxation parameter 
        // is not causing problems with time consistency
        UEqn = (ddtUEqn + HUEqn);
        UEqn.source() = S0;
        UEqn.boundaryCoeffs() = B0;


        // --- PISO loop
        while (pimple.correct())
        {
            p_().storePrevIter();   
                 
            AU_() = UEqn.A();
            HU_() = UEqn.H();

            U_() = HU_()/AU_();

            phi_() =
            (
                (fvc::interpolate(HU_())/fvc::interpolate(AU_()))
               & mesh().Sf()
            );

//            phi_() += fvc::ddtPhiCorr((1.0/AU_())(), rho_, U_(), phi_());

//            phi += fvc::ddtPhiCorr
//                (
//                    (rho_/AU_())(),
//                    U_(),
//                    phi_()
//                );

            // Global flux continuity
//            adjustPhi(phiHbyA_(), U_(), p_());

            if (closedVolume_)
            {
                label intPatchID_ = 
                    mesh().boundaryMesh().findPatchID("interfaceShadow");

                forAll(phi_().boundaryField(), patchI)
                {
                    if
                    (
                        !phi_().boundaryField()[patchI].coupled()
                     && patchI != intPatchID_
//                     && isA<zeroGradientFvPatchVectorField>
//                        (
//                            U_().boundaryField()[patchI]
//                        )
                    )
                    {
                        phi_().boundaryField()[patchI] ==
                        (
                            U_().boundaryField()[patchI]
                            & mesh().Sf().boundaryField()[patchI]
                        );
                    }
                }
            }

            if (!closedVolume_)
            {
                forAll(phi_().boundaryField(), patchI)
                {
                    if
                    (
                        !phi_().boundaryField()[patchI].coupled()
                     && !p_().boundaryField()[patchI].fixesValue()
                    )
                    {
                        phi_().boundaryField()[patchI] ==
                        (
                            U_().boundaryField()[patchI]
                            & mesh().Sf().boundaryField()[patchI]
                        );
                    }
                }
            }

            if (closedVolume_)
            {
                label intPatchID_ = 
                    mesh().boundaryMesh().findPatchID("interfaceShadow");

                phi_().boundaryField()[intPatchID_] =
                    (
                        U_().boundaryField()[intPatchID_]
                      & mesh().Sf().boundaryField()[intPatchID_]
                    );

                scalarField weights =
                    mag(phi_().boundaryField()[intPatchID_]);

                if(mag(gSum(weights)) > VSMALL)
                {
                    weights /= gSum(weights);
                }

                phi_().boundaryField()[intPatchID_] -=
                    weights*gSum(phi_().boundaryField()[intPatchID_]);

                phi_().boundaryField()[intPatchID_] +=
                    p_().boundaryField()[intPatchID_].snGrad()
                   *mesh().magSf().boundaryField()[intPatchID_]
                   /AU_().boundaryField()[intPatchID_];
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
                    fvm::laplacian(1.0/fvc::interpolate(AU_()), p_()) 
                 == fvc::div(phi_())
                );

                label pRefCell = 0;
                scalar pRefValue = 0.0;
                bool pNeedRef = false;
                bool procHasRef = false;

                // Find reference cell
                if (closedVolume_)
                {
                    point refPointi(mesh().solutionDict().subDict("PIMPLE").lookup("pRefPoint"));
                    label refCelli = mesh().findCell(refPointi);
                    label hasRef = (refCelli >= 0 ? 1 : 0);
                    label sumHasRef = returnReduce<label>(hasRef, sumOp<label>());

                    if (sumHasRef != 1)
                    {
                        FatalError<< "Unable to set reference cell for field "
                                << p_().name()
                                << nl << "    Reference point pRefPoint"
                                << " found on " << sumHasRef << " domains (should be one)"
                                << nl << exit(FatalError);
                    }

                    if (hasRef)
                    {
                        pRefCell = refCelli;
                        procHasRef = true;
                    }

                    pRefValue =
                        readScalar(mesh().solutionDict().subDict("PIMPLE").lookup("pRefValue"));

//                    if (pNeedRef && procHasRef)
                    {
                        pEqn.source()[pRefCell] +=
                            pEqn.diag()[pRefCell]*pRefValue;

                        pEqn.diag()[pRefCell] +=
                            pEqn.diag()[pRefCell];
                    }
                }

//                if (closedVolume_)
//                {
//                    point refPoint(mesh().solutionDict().subDict("PIMPLE").lookup("pRefPoint"));
//                    pRefCell_ = mesh().findCell(refPoint);

//                    pEqn.setReference(pRefCell_, pRefValue_);
//                }

                pEqn.solve
                (
                    mesh().solutionDict()
                    .solver(p_().select(pimple.finalInnerIter()))
                );

                if (pimple.finalNonOrthogonalIter())
                {
                    phi_() -= pEqn.flux();
                }                            
            }

            //- Pressure relaxation except for last corrector
            if (!pimple.finalIter())
            {
                p_().relax();
            }

            // Update of pressure gradient
            gradp_() = fvc::grad(p_());

            // Momentum corrector
            U_() -= fvc::grad(p_())/AU_();

            U_().correctBoundaryConditions();

            // Update of velocity gradient
            gradU_() = fvc::grad(U_());       
        }

        {
            mrfZones_.translationalMRFs().correctBoundaryVelocity(U_(), phi_());

            myTimeIndex_ = mesh().time().timeIndex();
        }
    }
}

// ************************************************************************* //
