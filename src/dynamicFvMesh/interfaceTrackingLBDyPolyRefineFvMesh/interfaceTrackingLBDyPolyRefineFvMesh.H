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
    Foam::interfaceTrackingLBDyPolyRefineFvMesh

Description
    Combination of interfaceTrackingFvMesh and 
    loadBalanceDynamicPolyRefinementFvMesh

SourceFiles
    interfaceTrackingLBDyPolyRefineFvMesh.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved
    Holger Marschall, TU Darmstadt.  All rights reserved
    Mohammed Elwardi Fadeli (elwardifadeli@gmail.com)
    Constantin Habes, TU Darmstadt.

\*---------------------------------------------------------------------------*/

#ifndef interfaceTrackingLBDyPolyRefineFvMesh_H
#define interfaceTrackingLBDyPolyRefineFvMesh_H

#include "loadBalanceDynamicPolyRefinementFvMesh.H"
#include "interfaceTrackingFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class interfaceTrackingLBDyPolyRefineFvMesh Declaration
\*---------------------------------------------------------------------------*/

class interfaceTrackingLBDyPolyRefineFvMesh
:
    public loadBalanceDynamicPolyRefinementFvMesh,
    public interfaceTrackingFvMesh
{
    // Private data

        // Copy control

            //- Disallow default bitwise copy construct
            interfaceTrackingLBDyPolyRefineFvMesh(const interfaceTrackingLBDyPolyRefineFvMesh&);

            //- Disallow default bitwise assignment
            void operator=(const interfaceTrackingLBDyPolyRefineFvMesh&);


public:

    //- Runtime type information
    TypeName("interfaceTrackingLBDyPolyRefineFvMesh");


    // Constructors

        //- Construct from IOobject and optional name for the subdictionary
        //  containing the refinement (and possibly other controls). The second
        //  parameter is useful if we derive another dynamic mesh class from
        //  this one instead of topoChangerFvMesh base class. VV, 17/May/2018.
        interfaceTrackingLBDyPolyRefineFvMesh
        (
            const IOobject& io
        );


    //- Destructor
    virtual ~interfaceTrackingLBDyPolyRefineFvMesh();


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

        //- Return access to interfaceTracking dictionary
        inline const dictionary& interfaceTrackingDict() const
        {
            return interfaceTrackingFvMesh::motionDict();
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
