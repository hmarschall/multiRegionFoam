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

Class
    Foam::epRegionModel.C

SourceFiles
    epRegionModel.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceNoneDWModel_H
#define surfaceNoneDWModel_H

#include "surfaceDissolvedModel.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace surfaceDissolvedModels
{

/*---------------------------------------------------------------------------*\
                           Class dissolved Declaration
\*---------------------------------------------------------------------------*/

class surfaceNoneDWModel
:
    public surfaceDissolvedModel
{
        //- Dictionary
        const dictionary& dict_;

public:

    //- Type name
        TypeName("surfaceNone");

    // Constructors

        surfaceNoneDWModel
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~surfaceNoneDWModel();


    // Member Functions
        //- update
        virtual void update(const word& clName);

        //- update dmdt at the faceCells
        virtual void updatePatch(const word& patchName);

        //- correct
        virtual void correct();

        //- solve, if necessary
        virtual void solve();

    //- IO

        virtual bool read(const dictionary& dict);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace surfaceDissolvedModels
} // End namespace Foam

#endif

// ************************************************************************* //
