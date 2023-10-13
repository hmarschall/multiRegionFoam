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

#include "fsiInterface.H"
#include "surfaceInterpolationScheme.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionInterfaces
{
    defineTypeNameAndDebug(fsiInterface, 0);

    addToRunTimeSelectionTable
    (
        regionInterfaceType,
        fsiInterface,
        IOdictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionInterfaces::fsiInterface::fsiInterface
(
    const word& type,
    const dictionary& dict,
    const Time& runTime,
    const fvPatch& patchA,
    const fvPatch& patchB
)
:
    regionInterfaceType(type, dict, runTime, patchA, patchB),
    dict_(dict.subDict(type + "Coeffs"))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionInterfaces::fsiInterface::makeInterfaceToInterface() const
{
    if (interfaceToInterfacePtr_.valid())
    {
        FatalErrorIn
        (
            "void Foam::interfaceToInterfaceMapping::"
            "makeInterfaceToInterface() const"
        )   << "Mapping object already set!" << abort(FatalError);
    }

    // Lookup the type
    const word type = interfaceProperties().lookupOrDefault<word>
    (
        "interfaceTransferMethod", "GGI"
    );

    // Assume meshA as solid region
    const fvMesh* solidMesh = &meshA();
    const fvMesh* fluidMesh = &meshB();

    // Check if meshB is solid mesh
    if (meshB().foundObject<pointVectorField>("pointD"))
    {
        solidMesh = &meshB();
        fluidMesh = &meshA();
    }

    Info<< endl << "Moving " << solidMesh->name()
        << " in the deformed configuration" << endl;

    // Get mesh points of solid mesh in initial configuration
    pointField newSolidMeshPoints = solidMesh->allPoints();

    // Get point displacement field of solid region
    const pointVectorField& pointD =
        solidMesh->lookupObject<pointVectorField>("pointD");

    // Add current displacement to mesh points
    newSolidMeshPoints += pointD;

    // Move the mesh in the deformed configuration
    const_cast<fvMesh&>(*solidMesh).movePoints(newSolidMeshPoints);

    Info<< "Creating updated interpolator" << endl;

    // Create interface to interface mapping
    interfaceToInterfacePtr_ =
    (
        interfaceToInterfaceMapping::New
        (
            type,
            interfaceProperties().subDict(type + "Coeffs"),
            meshA().boundaryMesh()[patchAID()],
            meshB().boundaryMesh()[patchBID()],
            globalPatchA(),
            globalPatchB()
        )
    );

    Info<< "Moving " << solidMesh->name()
        << " back in the initial configuration" << endl;

    // Substract displacement from mesh points again
    newSolidMeshPoints -= pointD;

    // Move the mesh back in the initial configuration
    const_cast<fvMesh&>(*solidMesh).movePoints(newSolidMeshPoints);
}

void Foam::regionInterfaces::fsiInterface::correct()
{
    // do nothing
}

Foam::scalar Foam::regionInterfaces::fsiInterface::getMinDeltaT()
{
    return GREAT;
}
// ************************************************************************* //
