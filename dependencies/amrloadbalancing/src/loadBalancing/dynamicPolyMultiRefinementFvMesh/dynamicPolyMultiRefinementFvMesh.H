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
    Foam::dynamicPolyMultiRefinementFvMesh

Description
    Adaptive mesh refinement for isotropic refinement of arbitrary polyhedral
    cells in 3D and arbitrary prismatic cells in 2D with composite refinement
    criteria.
    Automatically switches between:
    - 3D refinement uses polyhedralRefinement engine
    - 2D refinement uses prismatic2DRefinement engine

    Supports multiple refinement selection criteria.

SourceFiles
    dynamicPolyMultiRefinementFvMesh.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved
    Mohammed Elwardi Fadeli (elwardifadeli@gmail.com)

Notes
    Generalisation and refactorisation of dynamicRefineMesh for polyhedral cells
    in 3D and prismatic cells in 2D.

\*---------------------------------------------------------------------------*/

#ifndef dynamicPolyMultiRefinementFvMesh_H
#define dynamicPolyMultiRefinementFvMesh_H

#include "Tuple2.H"
#include "topoChangerFvMesh.H"
#include "refinementDetails.H"
#include "refinementSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class dynamicPolyMultiRefinementFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicPolyMultiRefinementFvMesh
:
    public virtual topoChangerFvMesh
{
    // Private data

        //- Dictionary containing dynamic mesh coefficient controls
        dictionary refinementDict_;

        //- List of refinements
        List<Tuple2<refinementDetails, autoPtr<refinementSelection>>> refinements_;
        
        //- Refinement interval (shared by all refinement ops)
        label refinementInterval_;

        //- Unrefinement interval (shared by all refinement ops)
        label unrefinementInterval_;

        // Copy control

            //- Disallow default bitwise copy construct
            dynamicPolyMultiRefinementFvMesh(const dynamicPolyMultiRefinementFvMesh&);

            //- Disallow default bitwise assignment
            void operator=(const dynamicPolyMultiRefinementFvMesh&);

        // Private member functions

            //- map a list of cell IDs after refinements
            void mapCellList
            (
                labelList& lst,
                const mapPolyMesh& topoChangeMap
            ) const;

            void mapPointList
            (
                labelList& lst,
                const mapPolyMesh& topoChangeMap
            ) const;
public:

    //- Runtime type information
    TypeName("dynamicPolyMultiRefinementFvMesh");


    // Constructors

        //- Construct from IOobject and optional name for the subdictionary
        //  containing the refinement (and possibly other controls). The second
        //  parameter is useful if we derive another dynamic mesh class from
        //  this one instead of topoChangerFvMesh base class. VV, 17/May/2018.
        dynamicPolyMultiRefinementFvMesh
        (
            const IOobject& io,
            const word subDictName = "dynamicRefineFvMesh"
        );


    //- Destructor
    virtual ~dynamicPolyMultiRefinementFvMesh();


    // Member Functions

        //- Return access to dictionary
        inline const dictionary& refinementDict() const
        {
            return refinementDict_;
        }

        //- Is this the first update in the current time-step?
        bool firstUpdate() const;

        //- Update the mesh for topology change
        virtual bool update();

        //- Single refinement mesh update
        bool singleRefUpdate();

        //- Multi refinement mesh update
        bool multiRefUpdate();

        //- Unrefinement mesh update
        bool unrefUpdate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
