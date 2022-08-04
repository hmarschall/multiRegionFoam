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

Class
    Foam::refinementDetails

Description
    Refinement parameters reader

SourceFiles
    refinementDetails.C

Author
    Mohammed Elwardi Fadeli (elwardifadeli@gmail.com)


\*---------------------------------------------------------------------------*/

#ifndef refinementDetails_H
#define refinementDetails_H

#include "topoChangerFvMesh.H"
#include "refinementSelection.H"
#include "enhancedRefinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class refinementDetails Declaration
\*---------------------------------------------------------------------------*/

class refinementDetails
{
    // Private data

        //- Dictionary containing dynamic mesh coefficient controls
        dictionary refinementDict_;

        //- Separate refinement/unrefinement: off by default, meaning that
        //  refinement and unrefinement can be performed in the same step
        Switch separateUpdates_;

        //- Current time index (helper variable to skip multiple topo changes in
        //  a single time step)
        label curTimeIndex_;

        //- Max refinement level
        label maxRefinementLevel_;

        //- Switch whether to use edge based consistency on refinement
        Switch edgeBasedConsistency_;

        //- Number of buffer layers for refinement
        label nRefinementBufferLayers_;

        //- Number of buffer layers for unrefinement, controlling how far
        //  the unrefinement region needs to be from current refinement
        //  region
        label nUnrefinementBufferLayers_;


public:

    //- Runtime type information
    TypeName("refinementDetails");


    // Constructors

        refinementDetails();

        //- Construct from IOobject and optional name for the subdictionary
        //  containing the refinement (and possibly other controls). The second
        //  parameter is useful if we derive another dynamic mesh class from
        //  this one instead of topoChangerFvMesh base class. VV, 17/May/2018.
        refinementDetails
        (
            const Time& time,
            const IOobject& io,
            const objectRegistry& db,
            const dictionary& refinementDict
        );


    //- Destructor
    virtual ~refinementDetails();


    // Member Functions

        //- Helper function for reading the dictionary and updating the data
        void readDict();

        //- Return refinement name
        inline word name() const {
            return refinementDict().name().removeRepeated('"').replaceAll("::", "");
        }

        //- Return access to dictionary
        inline const dictionary& refinementDict() const
        {
            return refinementDict_;
        }

        //- Always turn off seperate updates
        inline Switch separateUpdates() const
        {
            return false;
        }

        //- Return current time index
        inline label curTimeIndex() const
        {
            return curTimeIndex_;
        }

        //- Set current time index
        void setCurTimeIndex(label newCurTimeIndex)
        {
            curTimeIndex_ = newCurTimeIndex;
        }

        //- Set refinement settings from members
        void setRefinementSettings(enhancedRefinement&) const;

        //- Access to max refinement level
        label& maxRefinementLevel() {
            return maxRefinementLevel_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
