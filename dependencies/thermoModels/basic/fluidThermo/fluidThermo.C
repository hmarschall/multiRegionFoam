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

#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidThermo, 0);
    defineRunTimeSelectionTable(fluidThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidThermo::fluidThermo(const fvMesh& mesh, const word& phaseName)
:
    basicThermo(mesh, phaseName)
{}



// Foam::fluidThermo::fluidThermo
// (
//     const fvMesh& mesh,
//     const dictionary& dict,
//     const word& phaseName
// )
// :
//     basicThermo(mesh, dict, phaseName)
// {}


// Foam::fluidThermo::fluidThermo
// (
//     const fvMesh& mesh,
//     const word& phaseName,
//     const word& dictionaryName
// )
// :
//     basicThermo(mesh, phaseName, dictionaryName)
// {}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidThermo> Foam::fluidThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<fluidThermo>(mesh, phaseName);
}


// Foam::autoPtr<Foam::fluidThermo> Foam::fluidThermo::New
// (
//     const fvMesh& mesh,
//     const word& phaseName,
//     const word& dictName
// )
// {
//     return basicThermo::New<fluidThermo>(mesh, phaseName, dictName);
// }


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidThermo::~fluidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fluidThermo::nu() const
{
    return mu()/rho();
}


Foam::tmp<Foam::scalarField> Foam::fluidThermo::nu(const label patchi) const
{
    return mu(patchi)/rho(patchi);
}

// ************************************************************************* //
