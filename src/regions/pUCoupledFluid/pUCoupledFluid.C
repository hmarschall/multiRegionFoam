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
#include "pUCoupledFluid.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(pUCoupledFluid, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        pUCoupledFluid,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::pUCoupledFluid::pUCoupledFluid
(
    const Time& runTime,
    const word& regionName
)
:
    regionType(runTime, regionName),
    
    regionName_(regionName),
    U_(nullptr),
    phi_(nullptr),
    p_(nullptr),
    presSource_(nullptr),
    Up_(nullptr),
    rAU_(nullptr),
    gradp_(nullptr), 
    gradU_(nullptr), 
    pcorrTypes_(),   
    pcorr_(nullptr),
    pMin_
    (
        "pMin", 
        dimPressure,
        0
    ),
    pMax_
    (
        "pMax",
        dimPressure,
        0
    ),
    UMax_
    (
        "UMax",
        dimVelocity,
        0
    ),
    smallU_
    (
        "smallU",
        dimVelocity,
        1e-10
    ),
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
        readScalar
        (
            mesh().solutionDict().subDict("blockSolver").lookup("pRefValue")
        )
    ),
    whichProcHasRef_(Pstream::nProcs(), 0),
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
    ),
    firstUpdate_(true)
{
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

    // look up block vector pressure field from object registry
    if (mesh().foundObject<volVector4Field>("Up"))
    {
        Info << nl << "Using already existing block vector pressure field in region "
             << mesh().name()
             << " for regionType "
             << this->name()
             << nl << endl;

        Up_.reset
        (
            const_cast<volVector4Field*>
            (
                &mesh().lookupObject<volVector4Field>("Up")
            )
        );
    }
    else // create new block vector pressure field
    {
        Up_.reset
        (
            new volVector4Field
            (
                IOobject
                (
                    "Up",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedVector4("zero", dimless, vector4::zero)
            )
        );
    }

    // look up rAU field from object registry
    if (mesh().foundObject<volScalarField>("rAU"))
    {
        Info << nl << "Using already existing rAU field in region "
             << mesh().name()
             << " for regionType "
             << this->name()
             << nl << endl;

        rAU_.reset
        (
            const_cast<volScalarField*>
            (
                &mesh().lookupObject<volScalarField>("rAU")
            )
        );
    }
    else // create new rAU field
    {
        rAU_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "rAU",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                (dimVolume*dimTime)/dimMass,
                zeroGradientFvPatchScalarField::typeName
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
    else // create new viscosity field
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

    //- read field bounds
    dictionary fieldBounds = mesh().solutionDict().subDict("fieldBounds");
    fieldBounds.lookup(p_().name()) >> pMin_.value() >> pMax_.value();
    fieldBounds.lookup(U_().name()) >> UMax_.value();

    gradp_().checkIn();
    gradU_().checkIn();

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

    mesh().schemesDict().setFluxRequired(p_().name());

    #include "setRefCell.H"
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::pUCoupledFluid::~pUCoupledFluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::pUCoupledFluid::correct()
{
    if(firstUpdate_)
    {
        firstUpdate_=false;
    }
    else
    {
        #include "correctPhi.H"
    }
}

void Foam::regionTypes::pUCoupledFluid::setRDeltaT()
{
    #include "CourantNo.H"
    #include "setDeltaT.H"
}


void Foam::regionTypes::pUCoupledFluid::setCoupledEqns()
{   
    // Store p field for outer correction loop
    p_().storePrevIter();

    // Up equation name
    word UpEqnName = Up_().name() + mesh().name() + "Eqn";
    // Initialize the Up block system
    fvVector4Matrices.set
    (
        UpEqnName,
        new fvBlockMatrix<vector4>(Up_())
    );

    fvBlockMatrix<vector4>* UpEqn = getCoupledEqn<fvBlockMatrix,vector4>(UpEqnName);

    if (myTimeIndex_ < mesh().time().timeIndex())
    {
        mrfZones_.translationalMRFs().correctMRF();
    }

    // Make the fluxes relative to the mesh motion
    //fvc::makeRelative(phi_(), U_());
    phi_() == (phi_() - fvc::meshPhi(rho_(), U_()));

    // Assemble and insert momentum equation
    #include "UEqn.H"

    // Assemble and insert pressure equation
    #include "pEqn.H"

    // Assemble and insert coupling terms
    #include "couplingTerms.H"

    myTimeIndex_ = mesh().time().timeIndex();
        
}


void Foam::regionTypes::pUCoupledFluid::postSolve()
{
    // Retrieve solution
    word UpEqnName = Up_().name() + mesh().name() + "Eqn";
    fvVector4Matrices[UpEqnName]->retrieveSolution(0, U_().internalField());
    fvVector4Matrices[UpEqnName]->retrieveSolution(3, p_().internalField());

    Info<< "Pressure max: " << gMax(p_()) << " min: " << gMin(p_()) << " mean: " << gAverage(p_())
    << nl << "Velocity max: " << gMax(U_()) << " min: " << gMax(U_()) << " mean: " << gAverage(U_())
    << endl;

    U_().correctBoundaryConditions();
    p_().correctBoundaryConditions();

    word pEqnName = p_().name() + mesh().name() + "Eqn";
    phi_() = (fvc::interpolate(U_()) & mesh().Sf()) 
        + fvScalarMatrices[pEqnName]->flux() 
        + presSource_();

    #include "boundPU.H"

    gradp_() = fvc::grad(p_());
    gradU_() = fvc::grad(U_()); 

    mrfZones_.translationalMRFs().correctBoundaryVelocity(U_(), phi_());

    p_().relax();

}


void Foam::regionTypes::pUCoupledFluid::solveRegion()
{
    // do nothing, add as required
}

// ************************************************************************* //
