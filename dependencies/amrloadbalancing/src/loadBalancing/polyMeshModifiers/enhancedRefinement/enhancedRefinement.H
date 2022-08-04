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
    Foam::enhancedRefinement

Description
    Extends the refinement class with a proper interface to act on refinement
    settings, fixes some bugs

SourceFiles
    enhancedRefinement.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.
    Mohammed Elwardi Fadeli

\*---------------------------------------------------------------------------*/

#ifndef enhancedRefinement_H
#define enhancedRefinement_H

#include "refinement.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class enhancedRefinement Declaration
\*---------------------------------------------------------------------------*/

class enhancedRefinement
:
    public refinement
{
public:

    //- Runtime type information
    TypeName("enhancedRefinement");


    // Constructors

        //- Construct from dictionary
        enhancedRefinement
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyTopoChanger& mme
        );


    //- Destructor
    virtual ~enhancedRefinement();


    // Member Functions

        // Access

            //- Return maxRefinementLevel
            inline label maxRefinementLevel() const
            {
                return maxRefinementLevel_;
            }

            //- Set maxRefinementLevel
            inline void setMaxRefinementLevel(const label& maxRefLevel)
            {
                maxRefinementLevel_ = maxRefLevel;
            }

            //- Return edgeBasedConsistency
            inline Switch edgeBasedConsistency() const
            {
                return edgeBasedConsistency_;
            }

            //- Set edgeBasedConsistency
            inline void setEdgeBasedConsistency(const bool& edgeConsistency)
            {
                edgeBasedConsistency_ = edgeConsistency;
            }

            //- Return nRefinementBufferLayers
            inline label nRefinementBufferLayers() const
            {
                return nRefinementBufferLayers_;
            }

            //- Set nRefinementBufferLayers
            inline void setNRefinementBufferLayers(const label& nRefLayers)
            {
                nRefinementBufferLayers_ = nRefLayers;
            }

            //- Return nRefinementBufferLayers
            inline label nUnrefinementBufferLayers() const
            {
                return nUnrefinementBufferLayers_;
            }

            //- Set nUnrefinementBufferLayers
            inline void setNUnrefinementBufferLayers(const label& nUnrefLayers)
            {
                nUnrefinementBufferLayers_ = nUnrefLayers;
            }

            //- Access to cellsToRefine
            inline labelList cellsToRefine()
            {
                return cellsToRefine_;
            }

            // Return const access to cell Level
            const labelList& cellLevel() const {
                return cellLevel_;
            }

            // Shadows buggy refinement::faceConsistentRefinement
            label faceConsistentRefinement
            (
                boolList& cellsToRefine
            ) const;

            // Shadows buggy refinement::faceConsistentUnrefinement
            label faceConsistentUnrefinement
            (
                boolList& cellsToUnrefine
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
