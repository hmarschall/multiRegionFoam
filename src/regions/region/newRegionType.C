/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::regionType> Foam::regionType::New
(
    const Time& runTime,
    const word& name,
    const word& modelType
)
{
    word regionType(modelType);

    Info << "Selecting regionType: "
        << regionType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(regionType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("regionType::New")
           << "Unknown region type "
           << regionType << endl << endl
           << "Valid region types are : " << endl
           << dictionaryConstructorTablePtr_->sortedToc()
           << exit(FatalError);
    }

    return cstrIter()(runTime, name);
}

// ************************************************************************* //
