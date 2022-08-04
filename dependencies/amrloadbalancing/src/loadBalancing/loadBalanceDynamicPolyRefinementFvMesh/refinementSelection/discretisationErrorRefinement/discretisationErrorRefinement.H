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
    Foam::discretisationErrorRefinement

Description
    Selection of refinement cells based on discretisation error

SourceFiles
    discretisationErrorRefinement.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.
    Matthias Niethammer, TU Darmstadt.  All rights reserved.
    Franjo Juretic, Creative Fields Ltd.  All rights reserved.
    Constantin Habes, TU Darmstadt.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef discretisationErrorRefinement_H
#define discretisationErrorRefinement_H

#include "refinementSelection.H"
#include "errorEstimator.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class discretisationErrorRefinement Declaration
\*---------------------------------------------------------------------------*/

class discretisationErrorRefinement
:
    public refinementSelection
{
    // Private data

        //- Name of the field to refine
        word fieldName_;

        //- maximal discretisation error for the field
        scalar maxError_;

        //- Whether to use cell-point-cell interpolation to smooth out the field
        //  before selection. Switched off by default
        Switch cellPointCellSmoothing_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        discretisationErrorRefinement(const discretisationErrorRefinement&);

        //- Disallow default bitwise assignment
        void operator=(const discretisationErrorRefinement&);

public:

    //- Runtime type information
    TypeName("discretisationErrorRefinement");


    // Constructors

        //- Construct from components
        discretisationErrorRefinement
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~discretisationErrorRefinement();


    // Member Functions

        // Selection of refinement/unrefinement candidates

            //- Return transferable list of cells to refine
            virtual Xfer<labelList> refinementCellCandidates() const;

            //- Return transferable list of split points to unrefine
            virtual Xfer<labelList> unrefinementPointCandidates() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
