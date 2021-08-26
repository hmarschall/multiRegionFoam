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

#include "simpleControl.H"


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
//            this->time().timeName(),
//            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    k_(transportProperties_.lookup("k")),
    cv_(transportProperties_.lookup("cv")),
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
    ),
    coupledScalarEqnsTable_()
{}


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


void Foam::regionTypes::conductTemperature::solveCoupledPartitioned()
{
    // do nothing, add as required
}

void Foam::regionTypes::conductTemperature::solveRegion()
{
    Info << nl << "Solving for temperature in " << regionName_ << endl;
    simpleControl simpleControlRegion(*this);

    while (simpleControlRegion.correctNonOrthogonal())
    {
        tmp<fvScalarMatrix> TEqn
        (
            fvm::ddt(rho_*cv_, T_)
         ==
            fvm::laplacian(k_, T_, "laplacian(k,T)")
        );

        TEqn->relax();

        TEqn->solve();
    }
}

// HashTable
// <
//     Foam::tmp<Foam::fvScalarMatrix>, word, word::hash
// >
// Foam::regionTypes::conductTemperature::coupledScalarEqnsTable() const
// {
//     return coupledScalarEqnsTable_;
// }

// tmp<fvScalarMatrix> 
// Foam::regionTypes::conductTemperature::coupledFvScalarMatrix() const
// {
//     return tmp<fvScalarMatrix>
//     (
//         new fvScalarMatrix
//         (
//             fvm::ddt(rho_*cv_, T_)
//          ==
//             fvm::laplacian(k_, T_, "laplacian(k,T)")
//         )
//     );
// }


// ************************************************************************* //
