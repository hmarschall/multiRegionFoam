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
    Foam::regionInterfaceProperties

Description
    Simple class to hold region information for coupled region simulations.

    Gives per physics ('fluid', 'solid') the names of the regions. There
    is no assumption on this level that one region should only have one
    set of physics.

SourceFiles
    regionInterfaceProperties.C

\*---------------------------------------------------------------------------*/

#ifndef regionInterfaceProperties_H
#define regionInterfaceProperties_H

#include "foamTime.H"
#include "Pair.H"
#include "Tuple2.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class regionInterfaceProperties Declaration
\*---------------------------------------------------------------------------*/


class regionInterfaceProperties
:
    public HashTable<List<Tuple2<Pair<Pair<word> >, wordList > > >
{

        //- Disallow default bitwise copy construct
        regionInterfaceProperties(const regionInterfaceProperties&);

        //- Disallow default bitwise assignment
        void operator=(const regionInterfaceProperties&);

public:

    // Constructors

        //- Construct from components
        regionInterfaceProperties
        (
            const Time& runTime,
            const word
        );


    //- Destructor
    ~regionInterfaceProperties();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
