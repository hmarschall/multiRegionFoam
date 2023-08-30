    /*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "surfaceNoneDWModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceDissolvedModels
{
    defineTypeNameAndDebug(surfaceNoneDWModel, 0);

    addToRunTimeSelectionTable
    (
        surfaceDissolvedModel,
        surfaceNoneDWModel,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceDissolvedModels::surfaceNoneDWModel::surfaceNoneDWModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    surfaceDissolvedModel(mesh),

    dict_(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceDissolvedModels::surfaceNoneDWModel::~surfaceNoneDWModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::surfaceDissolvedModels::surfaceNoneDWModel::solve()
{
    // nothing
}


void Foam::surfaceDissolvedModels::surfaceNoneDWModel::correct()
{
    // nothing
}


void Foam::surfaceDissolvedModels::surfaceNoneDWModel::update
(
    const word& clName
)
{
    // nothing
}

void Foam::surfaceDissolvedModels::surfaceNoneDWModel::updatePatch
(
    const word& patchName
)
{
    // nothing
}

bool Foam::surfaceDissolvedModels::surfaceNoneDWModel::read(const dictionary& dict)
{
    return dict.isDict(type() + "Coeffs");
}
// ************************************************************************* //
