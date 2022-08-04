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

#include "interfaceTrackingLBDyPolyRefineFvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(interfaceTrackingLBDyPolyRefineFvMesh, 0);

    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        interfaceTrackingLBDyPolyRefineFvMesh,
        IOobject
    );

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceTrackingLBDyPolyRefineFvMesh::interfaceTrackingLBDyPolyRefineFvMesh
(
    const IOobject& io
)
:
    topoChangerFvMesh(io),
    loadBalanceDynamicPolyRefinementFvMesh(io),
    interfaceTrackingFvMesh(io)
{
    // Write mesh and modifiers
    topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
    topoChanger_.write();
    write();
    Info << "making Mesh" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceTrackingLBDyPolyRefineFvMesh::~interfaceTrackingLBDyPolyRefineFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::interfaceTrackingLBDyPolyRefineFvMesh::firstUpdate() const
{
    return loadBalanceDynamicPolyRefinementFvMesh::firstUpdate();
}


bool Foam::interfaceTrackingLBDyPolyRefineFvMesh::update()
{
    //Part 1 - move mesh
    bool hasMoved = interfaceTrackingFvMesh::update();

    //Part 2 - refine mesh and load balance
    // bool hasChanged = loadBalanceDynamicPolyRefinementFvMesh::update();
    bool hasChanged = dynamicPolyMultiRefinementFvMesh::update();

    return (hasChanged || hasMoved);
}


// ************************************************************************* //
