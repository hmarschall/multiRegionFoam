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

#include "dimensionedScalarFwd.H"
#include "heatTransferInterface.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionInterfaces
{
    defineTypeNameAndDebug(heatTransferInterface, 0);

    addToRunTimeSelectionTable
    (
        regionInterfaceType,
        heatTransferInterface,
        IOdictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionInterfaces::heatTransferInterface::heatTransferInterface
(
    const word& type,
    const dictionary& dict,
    const Time& runTime,
    const fvPatch& patchA,
    const fvPatch& patchB
)
:
    regionInterfaceType(type, dict, runTime, patchA, patchB),

    dict_(dict),

    K_(dimensionedScalar(dict_.lookup("K")))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dimensionedScalar& Foam::regionInterfaces::heatTransferInterface::K() const
{
    return K_;
}

void Foam::regionInterfaces::heatTransferInterface::correct()
{
    // do nothing
}

Foam::scalar Foam::regionInterfaces::heatTransferInterface::getMinDeltaT()
{
    return GREAT;
}
// ************************************************************************* //
