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
    const fvMesh& mesh,
    const word& regionName
)
:
    regionType(mesh, regionName),

    mesh_(mesh),
    regionName_(regionName),

    U_
    (
        IOobject
        (
            "U",
            this->time().timeName(),
            *this,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        *this
    ),
    phi_
    (
		IOobject
		(
			"phi",
            this->time().timeName(),
            *this,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		linearInterpolate(U_) & (*this).Sf()    
    ),

    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            this->time().constant(),
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    k_(transportProperties_.lookup("k")),
    cp_(transportProperties_.lookup("cp")),
    rho_(transportProperties_.lookup("rho")),
    T_
    (
        IOobject
        (
            "T",
            this->time().timeName(),
            *this,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        *this
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::transportTemperature::~transportTemperature()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::transportTemperature::correct()
{
    // do nothing, add as required
}


void Foam::regionTypes::transportTemperature::setRDeltaT()
{
    // do nothing, add as required
}


void Foam::regionTypes::transportTemperature::setCoupledEqns()
{
    dimensionedScalar alpha = k_/(rho_*cp_);

    fvScalarMatrix TEqn =
    (
        fvm::ddt(T_)
      + fvm::div(phi_, T_)
     ==
        fvm::laplacian(alpha, T_)
    );

    fvScalarMatrices.set
    (
        T_.name() + this->name() + "Eqn",
        new fvScalarMatrix(TEqn)
    );
}

void Foam::regionTypes::transportTemperature::solveRegion()
{
    // do nothing, add as required
}

// ************************************************************************* //
