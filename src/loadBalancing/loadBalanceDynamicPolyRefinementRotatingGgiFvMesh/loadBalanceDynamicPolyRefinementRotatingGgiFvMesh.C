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

\*---------------------------------------------------------------------------*/

#include "loadBalanceDynamicPolyRefinementRotatingGgiFvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(loadBalanceDynamicPolyRefinementRotatingGgiFvMesh, 0);

addToRunTimeSelectionTable
(
    topoChangerFvMesh,
    loadBalanceDynamicPolyRefinementRotatingGgiFvMesh,
    IOobject
);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::loadBalanceDynamicPolyRefinementRotatingGgiFvMesh::loadBalanceDynamicPolyRefinementRotatingGgiFvMesh
(
    const IOobject& io,
    const word refinementSubDictName,
    const word loadBalanceSubDictName,
    const word rotatingGgiSubDictName
)
:
    topoChangerFvMesh(io),
    loadBalanceDynamicPolyRefinementFvMesh(io, refinementSubDictName, loadBalanceSubDictName),
    rotatingGgiFvMesh(io, rotatingGgiSubDictName)
{
    // Write mesh and modifiers
    topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
    topoChanger_.write();
    write();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::loadBalanceDynamicPolyRefinementRotatingGgiFvMesh::~loadBalanceDynamicPolyRefinementRotatingGgiFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::loadBalanceDynamicPolyRefinementRotatingGgiFvMesh::firstUpdate() const
{
    return loadBalanceDynamicPolyRefinementFvMesh::firstUpdate();
}


bool Foam::loadBalanceDynamicPolyRefinementRotatingGgiFvMesh::update()
{
    
    //Part 1 - refine mesh and load balance
    bool hasChanged = loadBalanceDynamicPolyRefinementFvMesh::update();

    //Part 2 - move mesh
    bool hasMoved = rotatingGgiFvMesh::update();

    if (hasMoved)
    {
        Info << "Mesh has moved" << endl;
    }


    return (hasChanged || hasMoved);
}


// ************************************************************************* //
