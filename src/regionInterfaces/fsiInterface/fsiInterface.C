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
#include "scalar.H"
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
    dict_(dict.subDict(type + "Coeffs")),
    totalForce_(vector::zero)
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

    // Check if meshB is solid mesh
    if (meshB().foundObject<pointVectorField>("pointD"))
    {
        solidMesh = &meshB();
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

    // Set moving switch of mesh to be false again
    const_cast<fvMesh&>(*solidMesh).moving(false);
}

void
Foam::regionInterfaces::fsiInterface::setTotalForce(const vector& totalForce) const
{
    totalForce_ = totalForce;
}

void Foam::regionInterfaces::fsiInterface::correct()
{
    // do nothing
}

Foam::scalar Foam::regionInterfaces::fsiInterface::getMinDeltaT()
{
    return GREAT;
}

Foam::scalarField
Foam::regionInterfaces::fsiInterface::rawIntRes(const word& fldName) const
{
    //- RESIDUAL BASED ON TOTAL INTERFACE POINTS POSITIONS

    // Assume meshA as solid region
    bool solidMeshIsMeshA = true;
    const fvMesh* solidMesh = &meshA();
    const fvMesh* fluidMesh = &meshB();
    const fvPatch* solidPatch = &patchA();
    const fvPatch* fluidPatch = &patchB();

    // Check if meshB is solid mesh instead
    if (meshB().foundObject<pointVectorField>("pointD"))
    {
        solidMeshIsMeshA = false;
        solidMesh = &meshB();
        fluidMesh = &meshA();
        solidPatch = &patchB();
        fluidPatch = &patchA();
    }

    vectorField pResidual(patchA().patch().nPoints());

    const vectorField solidPatchPointD
    (
        solidMesh->lookupObject<pointVectorField>("pointD").internalField(),
        solidPatch->patch().meshPoints()
    );

    vectorField solidPatchPoints
    (
        solidMesh->allPoints(),
        solidPatch->patch().meshPoints()
    );

    // Location of solid patch points are the patch points
    // + their displacement
    solidPatchPoints = solidPatchPoints + solidPatchPointD;

    const vectorField fluidPatchPoints
    (
        fluidMesh->allPoints(),
        fluidPatch->patch().meshPoints()
    );

    if (solidMeshIsMeshA)
    {
        // Interpolate the fluid points to the solid side
        tmp<vectorField> tfluidPatchPointsToSolid =
            interpolatePointsFromB<vector>(fluidPatchPoints);

        vectorField fluidPatchPointsToSolid = tfluidPatchPointsToSolid();

        // Calculate residual
        pResidual = solidPatchPoints - fluidPatchPointsToSolid;
    }
    else
    {
        // Interpolate the solid points to the fluid side
        tmp<vectorField> tsolidPatchPointsToFluid =
            interpolatePointsFromB<vector>(solidPatchPoints);

        vectorField solidPatchPointsToFluid = tsolidPatchPointsToFluid();

        //Calculate residual
        pResidual = solidPatchPointsToFluid - fluidPatchPoints;
    }

    const primitivePatchInterpolation pfi(patchA().patch());

    tmp<vectorField> residual = pfi.pointToFaceInterpolate(pResidual);

    const tmp<scalarField> tmpRawResidual = mag(residual);
    const scalarField rawResidual = tmpRawResidual();

    return rawResidual;
}

Foam::scalar
Foam::regionInterfaces::fsiInterface::normIntRes(const word& fldName) const
{
    //- RESIDUAL BASED ON TOTAL INTERFACE POINTS POSITIONS

    // Assume meshA as solid region
    bool solidMeshIsMeshA = true;
    const fvMesh* solidMesh = &meshA();
    const fvMesh* fluidMesh = &meshB();
    const fvPatch* solidPatch = &patchA();
    const fvPatch* fluidPatch = &patchB();

    // Check if meshB is solid mesh instead
    if (meshB().foundObject<pointVectorField>("pointD"))
    {
        solidMeshIsMeshA = false;
        solidMesh = &meshB();
        fluidMesh = &meshA();
        solidPatch = &patchB();
        fluidPatch = &patchA();
    }

    vectorField residual(patchA().patch().nPoints());

    scalar n = 0;

    const vectorField solidPatchPointD
    (
        solidMesh->lookupObject<pointVectorField>("pointD").internalField(),
        solidPatch->patch().meshPoints()
    );

    vectorField solidPatchPoints
    (
        solidMesh->allPoints(),
        solidPatch->patch().meshPoints()
    );

    // Location of solid patch points are the patch points
    // + their displacement
    solidPatchPoints = solidPatchPoints + solidPatchPointD;

    const vectorField fluidPatchPoints
    (
        fluidMesh->allPoints(),
        fluidPatch->patch().meshPoints()
    );

    if (solidMeshIsMeshA)
    {
        // Interpolate the fluid points to the solid side
        tmp<vectorField> tfluidPatchPointsToSolid =
            interpolatePointsFromB<vector>(fluidPatchPoints);

        vectorField fluidPatchPointsToSolid = tfluidPatchPointsToSolid();

        // Calculate residual
        residual = solidPatchPoints - fluidPatchPointsToSolid;

        n = max
            (
                max
                (
                    Foam::sqrt(gSum(magSqr(solidPatchPoints))),
                    Foam::sqrt(gSum(magSqr(fluidPatchPointsToSolid)))
                ),
                SMALL
            );
    }
    else
    {
        // Interpolate the solid points to the fluid side
        tmp<vectorField> tsolidPatchPointsToFluid =
            interpolatePointsFromB<vector>(solidPatchPoints);

        vectorField solidPatchPointsToFluid = tsolidPatchPointsToFluid();

        //Calculate residual
        residual = solidPatchPointsToFluid - fluidPatchPoints;

        n = max
            (
                max
                (
                    Foam::sqrt(gSum(magSqr(solidPatchPointsToFluid))),
                    Foam::sqrt(gSum(magSqr(fluidPatchPoints)))
                ),
                SMALL
            );

    }

    //Return normalised residual
    return
    (
        Foam::sqrt(gSum(magSqr(residual)))/n
    );
}

void Foam::regionInterfaces::fsiInterface::info() const
{
    Info<< "Total force on solid interface: "
        << totalForce_<< endl;
}
// ************************************************************************* //
