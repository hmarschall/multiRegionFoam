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

    usesPIMPLE(true);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::regionTypeList::usesPIMPLE(const bool warn) const
{
    bool a = false;
    forAll(*this, i)
    {
        a = a || this->operator[](i).usesPIMPLE();
    }

    if (warn && this->size() && !a)
    {
        Info<< "No PIMPLE active" << endl;
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
    wordList updated;

    label n = 0;

    forAll(*this, i)
    {
        // mesh update (one sweep before solving)
        // Note: multiple coupled regions require an
        // updated system meshes prior to solution
        // (see Peric)
        if (findIndex(updated, this->operator[](i).mesh().name()) == -1)
        {
            updated.setSize(updated.size()+1);

            updated[n] = this->operator[](i).mesh().name();

            this->operator[](i).update();

            n++;
        }

        // correct properties
        this->operator[](i).correct();
    }
}


Foam::scalar Foam::regionTypeList::getMinDeltaT()
{
    scalar minDeltaT = GREAT;
    forAll(*this, i)
    {
        minDeltaT = min(minDeltaT, this->operator[](i).getMinDeltaT());
    }

    return minDeltaT;
}


void Foam::regionTypeList::solveRegion()
{
    forAll(*this, i)
    {
        // Solve for region-specific physics
        this->operator[](i).solveRegion();
    }
}

void Foam::regionTypeList::solvePIMPLE()
{
    // We do not have a top-level mesh. Construct the fvSolution for
    // the runTime instead.
    fvSolution solutionDict(runTime_);

    const dictionary& pimple = solutionDict.subDict("PIMPLE");

    int nOuterCorr(readInt(pimple.lookup("nOuterCorrectors")));

    //- PIMPLE loop
    for (int oCorr=0; oCorr<nOuterCorr; oCorr++)
    {
        forAll(*this, i)
        {
            this->operator[](i).prePredictor();
        }

        forAll(*this, i)
        {
            this->operator[](i).momentumPredictor();
        }

        forAll(*this, i)
        {
            this->operator[](i).pressureCorrector();
        }
    }
}

void Foam::regionTypeList::meshMotionCorrector()
{
    forAll(*this, i)
    {
        this->operator[](i).meshMotionCorrector();
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

void Foam::regionTypeList::postSolvePIMPLE()
{
    forAll(*this, i)
    {
        if (this->operator[](i).usesPIMPLE())
        {
            this->operator[](i).postSolve();
        }
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
