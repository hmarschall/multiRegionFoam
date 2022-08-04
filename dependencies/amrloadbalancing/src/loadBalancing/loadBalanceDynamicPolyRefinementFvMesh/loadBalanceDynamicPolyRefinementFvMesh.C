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

#include "dimensionedScalarFwd.H"
#include "dynamicFvMesh.H"
#include "dynamicPolyMultiRefinementFvMesh.H"
#include "loadBalanceDynamicPolyRefinementFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "loadBalanceFvMesh.H"
#include "refinementSelection.H"
#include "prismatic2DRefinement.H"
#include "polyhedralRefinement.H"
#include "topoChangerFvMesh.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(loadBalanceDynamicPolyRefinementFvMesh, 0);

addToRunTimeSelectionTable
(
    topoChangerFvMesh,
    loadBalanceDynamicPolyRefinementFvMesh,
    IOobject
);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::loadBalanceDynamicPolyRefinementFvMesh::loadBalanceDynamicPolyRefinementFvMesh
(
    const IOobject& io,
    const word refinementSubDictName,
    const word loadBalanceSubDictName
)
:
    topoChangerFvMesh(io),
    dynamicPolyMultiRefinementFvMesh(io, refinementSubDictName),
    loadBalanceFvMesh(io, loadBalanceSubDictName)
{
    // Write mesh and modifiers
    topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
    topoChanger_.write();
    write();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::loadBalanceDynamicPolyRefinementFvMesh::~loadBalanceDynamicPolyRefinementFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::loadBalanceDynamicPolyRefinementFvMesh::firstUpdate() const
{
    return dynamicPolyMultiRefinementFvMesh::firstUpdate();
}


bool Foam::loadBalanceDynamicPolyRefinementFvMesh::update()
{
    //Part 1 - update mesh
    bool hasChanged = dynamicPolyMultiRefinementFvMesh::update();

    //Part 2 - dynamic load balance
    if (time().timeIndex() > 1 && hasChanged && Pstream::parRun())
    {
        return loadBalanceFvMesh::loadBalance(loadBalanceDict());
    }

    return hasChanged;
}


// ************************************************************************* //
