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

Namespace
    Foam::globalFCParameters

Description
    global parameters for fuel cell modelling (Source: Vetter, Schumacher 2019)

\*---------------------------------------------------------------------------*/

#ifndef globalFCParameters_H
#define globalFCParameters_H

#include "scalar.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace regionTypes
{
    const dimensionedScalar pRef_("pRef", dimensionSet(1, -1, -2, 0, 0, 0, 0), 1e5);
    const dimensionedScalar TRef_("TRef", dimensionSet(0, 0, 0, 1, 0, 0, 0), 353.15);
    const dimensionedScalar RGas_("RGas", dimensionSet(1, 2, -2, -1, -1, 0, 0), 8.314472);
    const dimensionedScalar FConst_("FConst", dimensionSet(0, 0, 1, 0, -1, 1, 0), 9.64853e4);
    const scalar pi_(M_PI);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
