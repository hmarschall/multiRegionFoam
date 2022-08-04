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
    Foam::loadBalanceDynamicPolyRefinementRotatingGgiFvMesh

Description
    Adaptive mesh refinement for isotropic refinement of arbitrary polyhedral
    cells in 3D and arbitrary prismatic cells in 2D.
    Automatically switches between:
    - 3D refinement uses enhancedPolyhedralRefinement engine
    - 2D refinement uses enhancedPrismatic2DRefinement engine
    
    With load balancing capabilities

SourceFiles
    loadBalanceDynamicPolyRefinementRotatingGgiFvMesh.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved
    Holger Marschall, TU Darmstadt.  All rights reserved
    Mohammed Elwardi Fadeli (elwardifadeli@gmail.com)

Notes
    Generalisation and refactorisation of dynamicRefineMesh for polyhedral cells
    in 3D and prismatic cells in 2D including load balancing and multi-criteria
    refinement.

\*---------------------------------------------------------------------------*/

#ifndef loadBalanceDynamicPolyRefinementRotatingGgiFvMesh_H
#define loadBalanceDynamicPolyRefinementRotatingGgiFvMesh_H

#include "loadBalanceDynamicPolyRefinementFvMesh.H"
#include "rotatingGgiFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class loadBalanceDynamicPolyRefinementRotatingGgiFvMesh Declaration
\*---------------------------------------------------------------------------*/

class loadBalanceDynamicPolyRefinementRotatingGgiFvMesh
:
    public loadBalanceDynamicPolyRefinementFvMesh,
    public rotatingGgiFvMesh
{
    // Private data

        // Copy control

            //- Disallow default bitwise copy construct
            loadBalanceDynamicPolyRefinementRotatingGgiFvMesh(const loadBalanceDynamicPolyRefinementRotatingGgiFvMesh&);

            //- Disallow default bitwise assignment
            void operator=(const loadBalanceDynamicPolyRefinementRotatingGgiFvMesh&);


public:

    //- Runtime type information
    TypeName("loadBalanceDynamicPolyRefinementRotatingGgiFvMesh");


    // Constructors

        //- Construct from IOobject and optional name for the subdictionary
        //  containing the refinement (and possibly other controls). The second
        //  parameter is useful if we derive another dynamic mesh class from
        //  this one instead of topoChangerFvMesh base class. VV, 17/May/2018.
        loadBalanceDynamicPolyRefinementRotatingGgiFvMesh
        (
            const IOobject& io,
            const word refinementSubDictName = "dynamicRefineFvMesh",
            const word loadBalanceSubDictName = "loadBalanceFvMesh",
            const word rotatingGgiSubDictName = "rotatingGgiFvMesh"
        );


    //- Destructor
    virtual ~loadBalanceDynamicPolyRefinementRotatingGgiFvMesh();


    // Member Functions

        //- Return access to refinement dictionary
        inline const dictionary& refinementDict() const
        {
            return loadBalanceDynamicPolyRefinementFvMesh::refinementDict();
        }

        //- Return access to loadBalancing dictionary
        inline const dictionary& loadBalanceDict() const
        {
            return loadBalanceDynamicPolyRefinementFvMesh::loadBalanceDict();
        }


        //- Is this the first update in the current time-step?
        bool firstUpdate() const;

        //- Update the mesh for topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
