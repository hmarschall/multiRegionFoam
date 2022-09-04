/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

#include "regionTypeList.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "regionCouplePolyPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypeList::regionTypeList
(
    const Time& runTime
)
:
    PtrList<regionType>(),
    dict_
    (
        IOobject
        (
            "multiRegionProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    runTime_(runTime),
    region_(runTime)
{
    reset(region_);

    active(true);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::regionTypeList::active(const bool warn) const
{
    bool a = false;
    forAll(*this, i)
    {
        a = a || this->operator[](i).active();
    }

    if (warn && this->size() && !a)
    {
        Info<< "No models active" << endl;
    }

    return a;
}


void Foam::regionTypeList::reset(const regionProperties& rp)
{
    wordList regionTypes;

    label j = 0;

    forAll(rp, regionI)
    {
        const wordList& regions = rp[regionI].second();

        forAll(regions, regionI)
        {
            regionTypes.setSize(regionTypes.size()+1);
            regionTypes[j] = regions[regionI];

            j++;
        }
    }

    this->setSize(regionTypes.size());

    label i = 0;

    forAll(rp, regionI)
    {
        word meshName = rp[regionI].first();
        wordList regionTypes = rp[regionI].second();

        forAll(regionTypes, regionI)
        {
                Info<< "Creating region "
                    << meshName
                    << ": "
                    << regionTypes[regionI] 
                    << endl;

                this->set
                (
                    i++,
                    regionType::New
                    (
                        runTime_,
                        meshName,
                        regionTypes[regionI]
                    )
                );
        }
    }

//    wordList regionNames;

//    label j = 0;

//    forAllConstIter(HashTable<wordList>, rp, iter)
//    {
//        const wordList& regions = iter();

//        forAll(regions, regionI)
//        {
//            if (findIndex(regionNames, regions[regionI]))
//            {
//                regionNames.setSize(regionNames.size()+1);
//                regionNames[j] = regions[regionI];
//            }

//            j++;
//        }
//    }

//    this->setSize(regionNames.size());

//    label i = 0;

//    forAllConstIter(HashTable<wordList>, rp, iter)
//    {
//        const word& modelType = iter.key();
//        const wordList& regions = iter();

//        if (regions.size())
//        {
//            forAll(regions, regionI)
//            {
//                Info << "Creating " << regions[regionI] << endl;

//                this->set
//                (
//                    i++,
//                    regionType::New
//                    (
//                        runTime_,
//                        regions[regionI],
//                        modelType
//                    )
//                );
//            }
//        }
//    }

    // attach patches of regionCouplePolyPatch type
    forAll(*this, i)
    {
        dynamicFvMesh& mesh = const_cast<dynamicFvMesh&>(this->operator[](i).mesh());

        {
            const polyPatchList& patches = mesh.boundaryMesh();

            forAll (patches, patchI)
            {
                if (isType<regionCouplePolyPatch>(patches[patchI]))
                {
                    const regionCouplePolyPatch& rcp =
                        refCast<const regionCouplePolyPatch>(patches[patchI]);

                    // Attach it here, if slave
                    if (!rcp.master()) rcp.attach();
                }
            }

            // Force recalculation of weights
            mesh.surfaceInterpolation::movePoints();
        }
    }
}


void Foam::regionTypeList::preSolve()
{
    forAll(*this, i)
    {
        // mesh update (one sweep before solving)
        // Note: multiple coupled regions require an
        // updated system meshes prior to solution
        // (see Peric)
        this->operator[](i).update();

        // correct properties
        this->operator[](i).correct();
    }
}


void Foam::regionTypeList::setRDeltaT()
{
    forAll(*this, i)
    {
        this->operator[](i).setRDeltaT();
    }
}


void Foam::regionTypeList::solveRegion()
{
    for (int j=0; j<5; j++)
    {
        forAll(*this, i)
        {
            // Solve for region-specific physics
            // This might require outer loops if
            // coupling is achieved only by mutual
            // boundary condition updates
//            for (int j=0; j<5; j++)
            {
                this->operator[](i).solveRegion();
            }
        }
    }
}

void Foam::regionTypeList::setCoupledEqns()
{
    forAll(*this, i)
    {
        this->operator[](i).setCoupledEqns();
    }
}

void Foam::regionTypeList::postSolve()
{
    forAll(*this, i)
    {
        this->operator[](i).postSolve();
    }
}

void Foam::regionTypeList::clear()
{
    forAll(*this, i)
    {
        this->operator[](i).clear();
    }
}

// ************************************************************************* //
