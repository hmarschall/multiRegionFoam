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
#include "conductTemperature.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(conductTemperature, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        conductTemperature,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::conductTemperature::conductTemperature
(
    const fvMesh& mesh,
    const word& regionName
)
:
    regionType(mesh, regionName),

    mesh_(mesh),
    regionName_(regionName),

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
    cv_
    (
        IOobject
        (
            "cv",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this,
        dimensionedScalar(transportProperties_.lookup("cv"))
    ),
    rho_
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this,
        dimensionedScalar(transportProperties_.lookup("rho"))
    ),

//    k_
//    (
//        IOobject
//        (
//            "k",
//            this->time().timeName(),
//            *this,
//            IOobject::READ_IF_PRESENT,
////            IOobject::MUST_READ,
//            IOobject::NO_WRITE
//        ),
//        *this,
//        dimensionedScalar(transportProperties_.lookup("k"))
//    ),
//    T_
//    (
//        IOobject
//        (
//            "T",
//            this->time().timeName(),
//            *this,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        *this,
//        dimensionedScalar("T0", dimTemperature, pTraits<scalar>::zero),
//        zeroGradientFvPatchScalarField::typeName
//    ),
    k_(nullptr),
    kSolid_(nullptr),
    T_(nullptr)
{
    k_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "k",
                this->time().timeName(),
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            *this,
            dimensionedScalar("k", dimensionSet(0,2,-1,0,0,0,0), 0)
        )
    );

    kSolid_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "kSolid",
                this->time().timeName(),
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            *this,
            dimensionedScalar(transportProperties_.lookup("kSolid"))
        )
    );

    T_.reset
    (
        new volScalarField
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
    );

    k_() = kSolid_()/(rho_*cv_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::conductTemperature::~conductTemperature()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::conductTemperature::correct()
{
    // do nothing, add as required
}


void Foam::regionTypes::conductTemperature::setRDeltaT()
{
    // do nothing, add as required
}


void Foam::regionTypes::conductTemperature::setCoupledEqns()
{
    fvScalarMatrix TEqn =
    (
        fvm::ddt(rho_*cv_, T_())
     ==
        fvm::laplacian(kSolid_(), T_(), "laplacian(k,T)")
    );

    fvScalarMatrices.set
    (
        T_().name() + this->name() + "Eqn",
        new fvScalarMatrix(TEqn)
    );
}

void Foam::regionTypes::conductTemperature::solveRegion()
{
    // do nothing, add as required
}


// ************************************************************************* //
