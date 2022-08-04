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
    Foam::regionsRefinement

Description
    Selection of refinement cells contained with any cellSource region.
    unrefinementCandidates is set to zero since the mesh remains static.
    Set dynamicFvMesh to staticFvMesh during simulation run for better
    performance.

Usage
    Example usage in dynamicMeshDict.refinements:
    Any cellSource type can be define, e.g. boxToCell, sphereToCell, etc.
    \verbatim
        refinementSelection
        {
            type        regionsRefinement;
            regions
            (
                boxToCell
                {
                    minLevel 3;
                    box (Xmin Ymin Zmin)(Xmax Ymax Zmax);
                }
            );
              cellPointCellSmoothing on;
            }
    \endverbatim

SourceFiles
    regionsRefinement.C

Author
    Daniel Deising, TU Darmstadt. All rights reserved.
    Milad Bagheri, TU Darmstadt. All rights reserved.
    
\*---------------------------------------------------------------------------*/

#ifndef regionsRefinement_H
#define regionsRefinement_H

#include "refinementSelection.H"
#include "volFields.H"
#include "Switch.H"
#include "topoSetSource.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class regionsRefinement Declaration
\*---------------------------------------------------------------------------*/

class regionsRefinement
:
    public refinementSelection
{
    // Private data

        //- Reference to fvMesh
        const fvMesh& mesh_;
        
        //- Regions can be set to be quasi-static refined
        PtrList<entry> refinedRegions_;

        //- Whether to use cell-point-cell interpolation to smooth
        //  out the field before selection. Switched off by default
        Switch cellPointCellSmoothing_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionsRefinement(const regionsRefinement&);

        //- Disallow default bitwise assignment
        void operator=(const regionsRefinement&);

public:

    //- Runtime type information
    TypeName("regionsRefinement");


    // Constructors

        //- Construct from components
        regionsRefinement(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~regionsRefinement();


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
