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
#include "transportTemperature.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(transportTemperature, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        transportTemperature,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::transportTemperature::transportTemperature
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
    k_(transportProperties_.lookup("k")),
    cp_(transportProperties_.lookup("cp")),
    rho_(transportProperties_.lookup("rho")),

    U_(nullptr),
    alpha_(nullptr),
    phi_(nullptr),
    T_(nullptr)
{
    // set velocity field
    // Postponing field creation since U is probably provided by
    // another regionType, e.g. icoFluid, and thus to be re-used.
    U_ = lookupOrRead<volVectorField>(mesh(), "U");

    // set flux field
    phi_ = lookupOrRead<surfaceScalarField>
    (
        mesh(),
        "phi",
        true,
        false,
        linearInterpolate(U_()) & mesh().Sf()
    );

    // set thermal diffusivity field
    alpha_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "alpha", 
        k_/(rho_*cp_),
        false
    );

    // set temperature field
    T_ = lookupOrRead<volScalarField>(mesh(), "T");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::transportTemperature::~transportTemperature()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::transportTemperature::correct()
{
    // do nothing, add as required
}


Foam::scalar Foam::regionTypes::transportTemperature::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::transportTemperature::setCoupledEqns()
{
    fvScalarMatrix TEqn =
    (
        fvm::ddt(T_())
      + fvm::div(phi_(), T_())
     ==
        fvm::laplacian(alpha_(), T_())
    );

    fvScalarMatrices.set
    (
        T_().name() + mesh().name() + "Eqn",
        new fvScalarMatrix(TEqn)
    );
}

void Foam::regionTypes::transportTemperature::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportTemperature::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportTemperature::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportTemperature::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportTemperature::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportTemperature::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //
