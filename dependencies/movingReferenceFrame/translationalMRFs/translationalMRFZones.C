/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "translationalMRFZones.H"
#include "foamTime.H"
#include "fvMesh.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<translationalMRFZone>, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::translationalMRFZones::translationalMRFZones(const fvMesh& mesh)
:
    IOPtrList<translationalMRFZone>
    (
        IOobject
        (
            "translationalMRFZones",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        translationalMRFZone::iNew(mesh)
    ),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::translationalMRFZones::correctMRF()
{
    forAll (*this, i)
    {
        operator[](i).correctMRF();
    }
}


void Foam::translationalMRFZones::addFrameAcceleration
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho
)
{
    forAll (*this, i)
    {
        operator[](i).addFrameAcceleration(UEqn, rho);
    }
}

void Foam::translationalMRFZones::addFrameAcceleration
(
    fvVectorMatrix& UEqn
)
{
    forAll (*this, i)
    {
        operator[](i).addFrameAcceleration(UEqn);
    }
}


void Foam::translationalMRFZones::correctBoundaryVelocity
(
    volVectorField& U,
    surfaceScalarField& phi
)
{
    forAll (*this, i)
    {
        operator[](i).correctBoundaryVelocity(U, phi);
    }
}

void Foam::translationalMRFZones::writeRestart()
{
    forAll (*this, i)
    {
        operator[](i).writeRestart();
    }
}

// ************************************************************************* //
