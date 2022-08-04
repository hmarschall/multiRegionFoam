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
    nuFluid_
    (
        transportProperties_.subDict(regionName_).lookup("nu")
    ),
    nu_(nullptr),
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
    adjustTimeStep_
    (
        mesh().time().controlDict()
        .lookupOrDefault<Switch>("adjustTimeStep", false)
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
    )
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
                mesh().time().deltaT()
            )
        );
    }

    // look up viscosity field from object registry
    if (mesh().foundObject<volScalarField>("nu"))
    {
        Info << nl << "Using already existing viscosity field in region "
             << mesh().name()
             << " for regionType "
             << this->name()
             << nl << endl;

        nu_.reset
        (
            const_cast<volScalarField*>
            (
                &mesh().lookupObject<volScalarField>("nu")
            )
        );
    }
    else // create new viscosity field
    {
        nu_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "nu",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                nuFluid_
            )
        );
    }

    //- read field bounds
    dictionary fieldBounds = mesh().solutionDict().subDict("fieldBounds");
    fieldBounds.lookup(p_().name()) >> pMin_.value() >> pMax_.value();
    fieldBounds.lookup(U_().name()) >> UMax_.value();

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

    // const objectRegistry& dbParent = mesh().thisDb().parent();
    // const objectRegistry& db = mesh().thisDb();

    // Info << nl <<"Objects registered to parent of region: " << mesh().name() << nl << endl;

    // forAllConstIter(HashTable<regIOobject*>, dbParent, iter)
    // {
    //     Info << " name : " << iter()->name() << nl
    //         << " type : " << iter()->type() << nl
    //         << endl;
    // }

    // Info << nl <<"Objects registered to mesh of region: " << mesh().name() << nl << endl;

    // forAllConstIter(HashTable<regIOobject*>, db, iter)
    // {
    //     Info << " name : " << iter()->name() << nl
    //         << " type : " << iter()->type() << nl
    //         << endl;
    // }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::pUCoupledFluid::~pUCoupledFluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::pUCoupledFluid::correct()
{
    if (mesh().changing())
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

    fvBlockMatrix<vector4>& UpEqn = getCoupledEqn<fvBlockMatrix,vector4>(UpEqnName);

    if (myTimeIndex_ < mesh().time().timeIndex())
    {
        mrfZones_.translationalMRFs().correctMRF();
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi_(), U_());

    // Assemble and insert momentum equation
    #include "UEqn.H"

    // Assemble and insert pressure equation
    #include "pEqn.H"

    // Assemble and insert coupling terms
    #include "couplingTerms.H"

    myTimeIndex_ = mesh().time().timeIndex();
        
}


void Foam::regionTypes::pUCoupledFluid::updateFields()
{
    // Retrieve solution
    word UpEqnName = Up_().name() + mesh().name() + "Eqn";
    fvVector4Matrices[UpEqnName]->retrieveSolution(0, U_().internalField());
    fvVector4Matrices[UpEqnName]->retrieveSolution(3, p_().internalField());

    U_().correctBoundaryConditions();
    p_().correctBoundaryConditions();

    word pEqnName = p_().name() + mesh().name() + "Eqn";
    phi_() = (fvc::interpolate(U_()) & mesh().Sf()) 
        + fvScalarMatrices[pEqnName]->flux() 
        + presSource_();

    #include "boundPU.H"

    mrfZones_.translationalMRFs().correctBoundaryVelocity(U_(), phi_());

    p_().relax();

}


void Foam::regionTypes::pUCoupledFluid::solveRegion()
{
    // do nothing, add as required
}

// ************************************************************************* //
