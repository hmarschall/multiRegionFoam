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
#include "scalar.H"

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
            false
        )
    ),
    meshPtr_(nullptr)
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

    if (mesh().solutionDict().found("PICARD"))
    {
        const dictionary& picardControls = mesh().solutionDict().subDict("PICARD");

        forAllConstIter(dictionary, picardControls, iter)
        {
            if (iter().isDict())
            {
                const dictionary& subFieldDict(iter().dict());
                word fldName= iter().keyword();

                maxCorr_.insert
                (
                    fldName,
                    readScalar(subFieldDict.lookup("maxCorr"))
                );

                relativeTolerance_.insert
                (
                    fldName,
                    readScalar(subFieldDict.lookup("relTol"))
                );

                convergenceTolerance_.insert
                (
                    fldName,
                    readScalar(subFieldDict.lookup("tolerance"))
                );
            }
        }
    }
}

Foam::scalar Foam::regionType::maxCorr(word name)
{
    HashTable<scalar>::iterator it = maxCorr_.find(name);

    if (it == maxCorr_.end()) // not found
    {
        // if
        // (
        //     relativeTolerance(name) != VGREAT
        //  || convergenceTolerance(name)!= VGREAT)
        // {
        //     FatalErrorIn("scalar Foam::regionType::maxCorr(word name)")
        //         << "Either relTol or tolerance are specified for the solution of field "
        //         << name << " in region " << mesh().name()
        //         << "but no number for maxCorr is provided"
        //         << abort(FatalError);
        // }

        // If not specified no picard iterations are applied
        return 0;
    }

    return *it;
}


Foam::scalar Foam::regionType::relativeTolerance(word name)
{
    HashTable<scalar>::iterator it = relativeTolerance_.find(name);

    if (it == relativeTolerance_.end()) // not found
    {

        // if
        // (
        //     maxCorr(name) != 0
        // )
        // {
        //     Warning
        //         << "A number for maxCorr is specified for the solution of field "
        //         << name << " in region " << mesh().name()
        //         << "but no number for relTol is provided." << nl
        //         << "Setting relTol to 0 by default" << nl
        //         << endl;

        //     return 0;
        // }

        // If not specified set tolerance to be satisfied by default
        return VGREAT;
    }

    return *it;
}

Foam::scalar Foam::regionType::convergenceTolerance(word name)
{
    HashTable<scalar>::iterator it = convergenceTolerance_.find(name);

    if (it == convergenceTolerance_.end()) // not found
    {
        // if
        // (
        //     maxCorr(name) != 0
        // )
        // {
        //     Warning
        //         << "A number for maxCorr is specified for the solution of field "
        //         << name << " in region " << mesh().name()
        //         << "but no number for tolerance is provided." << nl
        //         << "Setting tolerance to 0 by default" << nl
        //         << endl;

        //     return 0;
        // }

        // If not specified set tolerance to be satisfied by default
        return VGREAT;
    }

    return *it;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "regionTypeTemplates.C"
#endif

// ************************************************************************* //
