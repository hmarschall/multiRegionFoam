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

#include "regionType.H"
#include "multiRegionSystem.H"
#include "IOReferencer.H"

namespace Foam
{
    defineTypeNameAndDebug(regionType, 0);
    defineRunTimeSelectionTable(regionType, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::regionType::regionType
(
    const Time& runTime,
    const word& regionName
)
:
    IOdictionary
    (
        IOobject
        (
            regionName + "Dict",
//            // If region == "region0" then read from the main case
//            // Otherwise, read from the region/sub-mesh directory
//            bool(regionName == dynamicFvMesh::defaultRegion)
//          ? fileName(runTime.caseConstant())
//          : fileName(runTime.caseConstant()/regionName),
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true //true        // Needs to be saved in registry, since I am searching for regionTypes in Electrochem!
        )
    ),
    meshPtr_(nullptr)
    // meshPtr_
    // (
    //     dynamicFvMesh::New
    //     (
    //         IOobject
    //         (
    //             regionName,
    //             runTime.timeName(),
    //             runTime,
    //             IOobject::MUST_READ
    //         )
    //     )
    // ),
    // faceRegionAddressingIO_
    // (
    //     IOobject
    //     (
    //         "faceRegionAddressing",
    //         meshPtr_->time().findInstance(meshPtr_->meshDir(), "faces"),
    //         polyMesh::meshSubDir,
    //         meshPtr_(),
    //         IOobject::MUST_READ
    //     )
    // ),

    // cellMapIO_
    // (
    //     IOobject
    //     (
    //         "cellRegionAddressing",
    //         meshPtr_->time().findInstance(meshPtr_->meshDir(), "faces"),
    //         polyMesh::meshSubDir,
    //         meshPtr_(),
    //         IOobject::MUST_READ
    //     )
    // ),
    // patchesMapIO_
    // (
    //     IOobject
    //     (
    //         "boundaryRegionAddressing",
    //         meshPtr_->time().findInstance(meshPtr_->meshDir(), "faces"),
    //         polyMesh::meshSubDir,
    //         meshPtr_(),
    //         IOobject::MUST_READ
    //     )
    // ),
    // faceMap_(faceRegionAddressingIO_.size(), 1),
    // faceMask_(faceRegionAddressingIO_.size(), 1)
{
    // look up mesh from object registry
    if (runTime.foundObject<dynamicFvMesh>(regionName))
    {
        meshPtr_.reset
        (
            const_cast<dynamicFvMesh*>
            (
                &runTime.lookupObject<dynamicFvMesh>(regionName)
            )
        );
    }
    // or create new mesh
    else
    {
        meshPtr_ = dynamicFvMesh::New
        (
            IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        );
    }

    // forAll(faceRegionAddressingIO_, i)
    // {
    //     faceRegionAddressing_.insert
    //     (
    //         faceRegionAddressingIO_[i],
    //         i
    //     );
    // }

    // forAll(cellMapIO_, i)
    // {
    //     cellMap_.insert
    //     (
    //         cellMapIO_[i],
    //         i
    //     );
    // }

    // forAll(patchesMapIO_, i)
    // {
    //     patchesMap_.insert
    //     (
    //         patchesMapIO_[i],
    //         i
    //     );
    // }

    // forAll(faceMap_, i)
    // {
    //     faceMap_[i] = mag(faceRegionAddressingIO_[i]) - 1;
    //     faceMask_[i] = sign(faceRegionAddressingIO_[i]);
    // }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "regionTypeTemplates.C"
#endif

// ************************************************************************* //
