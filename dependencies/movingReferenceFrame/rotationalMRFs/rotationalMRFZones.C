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

#include "rotationalMRFZones.H"
#include "foamTime.H"
#include "fvMesh.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<rotationalMRFZone>, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rotationalMRFZones::rotationalMRFZones(const fvMesh& mesh)
:
    IOPtrList<rotationalMRFZone>
    (
        IOobject
        (
            "rotationalMRFZones",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        rotationalMRFZone::iNew(mesh)
    ),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::rotationalMRFZones::omega() const
{
    tmp<volVectorField> trotationalMRFZonesOmega
    (
        new volVectorField
        (
            IOobject
            (
                "rotationalMRFZonesOmega",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimless/dimTime, vector::zero)
        )
    );
    volVectorField& rotationalMRFZonesOmega = trotationalMRFZonesOmega();

    forAll (*this, i)
    {
        operator[](i).addOmega(rotationalMRFZonesOmega);
    }

    return trotationalMRFZonesOmega;
}


Foam::tmp<Foam::surfaceScalarField> Foam::rotationalMRFZones::fluxCorrection() const
{
    tmp<surfaceScalarField> trotationalMRFZonesPhiCorr
    (
        new surfaceScalarField
        (
            IOobject
            (
                "rotationalMRFZonesPhiCorr",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimVelocity*dimArea, 0)
        )
    );
    surfaceScalarField& rotationalMRFZonesPhiCorr = trotationalMRFZonesPhiCorr();

    forAll (*this, i)
    {
        operator[](i).relativeFlux(rotationalMRFZonesPhiCorr);
    }

    return trotationalMRFZonesPhiCorr;
}


Foam::tmp<Foam::surfaceScalarField> Foam::rotationalMRFZones::meshPhi() const
{
    tmp<surfaceScalarField> trotationalMRFZonesFaceU
    (
        new surfaceScalarField
        (
            IOobject
            (
                "rotationalMRFZonesFaceU",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimVolume/dimTime, 0)
        )
    );
    surfaceScalarField& rotationalMRFZonesFaceU = trotationalMRFZonesFaceU();

    forAll (*this, i)
    {
        operator[](i).meshPhi(rotationalMRFZonesFaceU);
    }

    return trotationalMRFZonesFaceU;
}


void Foam::rotationalMRFZones::addCoriolis(fvVectorMatrix& UEqn) const
{
    forAll (*this, i)
    {
        operator[](i).addCoriolis(UEqn);
    }
}


void Foam::rotationalMRFZones::relativeFlux(surfaceScalarField& phi) const
{
    forAll (*this, i)
    {
        operator[](i).relativeFlux(phi);
    }
}


void Foam::rotationalMRFZones::absoluteFlux(surfaceScalarField& phi) const
{
    forAll (*this, i)
    {
        operator[](i).absoluteFlux(phi);
    }
}


void Foam::rotationalMRFZones::addCoriolis
(
    const volScalarField& rho,
    fvVectorMatrix& UEqn
) const
{
    forAll (*this, i)
    {
        operator[](i).addCoriolis(rho, UEqn);
    }
}


void Foam::rotationalMRFZones::relativeFlux
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    forAll (*this, i)
    {
        operator[](i).relativeFlux(rho, phi);
    }
}


void Foam::rotationalMRFZones::absoluteFlux
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    forAll (*this, i)
    {
        operator[](i).absoluteFlux(rho, phi);
    }
}


void Foam::rotationalMRFZones::relativeVelocity(volVectorField& U) const
{
    forAll (*this, i)
    {
        operator[](i).relativeVelocity(U);
    }
}


void Foam::rotationalMRFZones::absoluteVelocity(volVectorField& U) const
{
    forAll (*this, i)
    {
        operator[](i).absoluteVelocity(U);
    }
}


void Foam::rotationalMRFZones::correctBoundaryVelocity(volVectorField& U) const
{
    forAll (*this, i)
    {
        operator[](i).correctBoundaryVelocity(U);
    }
}


Foam::tmp<Foam::volScalarField> Foam::rotationalMRFZones::Su
(
    const volScalarField& phi
) const
{
    tmp<volScalarField> tPhiSource
    (
        new volScalarField
        (
            IOobject
            (
                phi.name() + "Source",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", phi.dimensions()/dimTime, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& source = tPhiSource();

    // Due to gradient cacheing, must take a tmp field
    // HJ, 22/Apr/2016
    tmp<volVectorField> tgradPhi = fvc::grad(phi);
    const volVectorField& gradPhi = tgradPhi();

    forAll (*this, i)
    {
        operator[](i).Su(phi, gradPhi, source);
    }

    return tPhiSource;
}


Foam::tmp<Foam::volVectorField> Foam::rotationalMRFZones::Su
(
    const volVectorField& phi
) const
{
    tmp<volVectorField> tPhiSource
    (
        new volVectorField
        (
            IOobject
            (
                phi.name() + "Source",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", phi.dimensions()/dimTime, vector::zero),
            zeroGradientFvPatchVectorField::typeName
        )
    );
    volVectorField& source = tPhiSource();

    // Due to gradient cacheing, must take a tmp field
    // HJ, 22/Apr/2016
    tmp<volTensorField> tgradPhi = fvc::grad(phi);
    const volTensorField& gradPhi = tgradPhi();

    forAll (*this, i)
    {
        operator[](i).Su(phi, gradPhi, source);
    }

    return tPhiSource;
}


// ************************************************************************* //
