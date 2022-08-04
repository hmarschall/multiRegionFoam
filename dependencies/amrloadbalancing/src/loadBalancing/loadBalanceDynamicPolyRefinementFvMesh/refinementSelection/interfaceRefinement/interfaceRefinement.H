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
    Foam::interfaceRefinement

Description
    Selection of refinement cells based on a phase interface

SourceFiles
    interfaceRefinement.C

Author

\*---------------------------------------------------------------------------*/

#ifndef interfaceRefinement_H
#define interfaceRefinement_H

#include "refinementSelection.H"
#include "volFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class interfaceRefinement Declaration
\*---------------------------------------------------------------------------*/

class interfaceRefinement
:
    public refinementSelection
{
    // Private data

        //- Field to refine
        wordList fldNames_;

        //- Reference to fvMesh
        const fvMesh& mesh_;

        //- Number if refinement layers inside of the interface
        scalar innerRefLayers_;

        //- Number if refinement layers outside of the interface
        scalar outerRefLayers_;

        //- Whether to use cell-point-cell interpolation to smooth
        //  out the field before selection. Switched off by default
        Switch cellPointCellSmoothing_;

        //- Tolerance value to select interface cells
        // if mag(alpha_N - alpha_P) > tolerance_, 
        // cells N,P are on the interface
        scalar tolerance_;

        //- list of phase indicators
        mutable PtrList<volScalarField> isInterface_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        interfaceRefinement(const interfaceRefinement&);

        //- Disallow default bitwise assignment
        void operator=(const interfaceRefinement&);

        //- Update of interface refinement candidates
        virtual void interfaceCandidates() const;

public:

    //- Runtime type information
    TypeName("interfaceRefinement");


    // Constructors

        //- Construct from components
        interfaceRefinement(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~interfaceRefinement();


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
