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
    Foam::loadBalanceFvMesh

Description
    Load balancing fvMesh

SourceFiles
    loadBalanceFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef loadBalanceFvMesh_H
#define loadBalanceFvMesh_H

#include "topoChangerFvMesh.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class loadBalanceFvMesh Declaration
\*---------------------------------------------------------------------------*/

class loadBalanceFvMesh
:
    virtual public topoChangerFvMesh
{
protected:

    // Protected data

        //- Dictionary of motion control parameters
        dictionary dict_;

        //- Imbalance trigger 0 < imbalance < 1
        scalar imbalanceTrigger_;

    // Protected member functions

        // Helper functions 

            //- Reading the dictionary and updating the data
            void readDict();

    // Private Member Functions

        //- Disallow default bitwise copy construct
        loadBalanceFvMesh(const loadBalanceFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const loadBalanceFvMesh&);

        //- Check balance.  Return true if re-balancing is required
        bool checkLoadBalance(const scalarField& weights) const;


public:

    //- Runtime type information
    TypeName("loadBalanceFvMesh");


    // Constructors

        //- Construct from IOobject
        loadBalanceFvMesh
        (
            const IOobject& io,
            const word& subDictName = typeName
        );


    //- Destructor
    virtual ~loadBalanceFvMesh();


    // Member Functions

        //- Return access to dictionary
        inline const dictionary& dict() const
        {
            return dict_;
        }

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
