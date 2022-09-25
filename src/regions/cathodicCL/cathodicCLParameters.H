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
    Foam::catalystLayerParameters

Description
    catalyst layer parameters for fuel cell modelling (Source: Vetter, Schumacher 2019)

\*---------------------------------------------------------------------------*/

#ifndef catalystLayerParameters_H
#define catalystLayerParameters_H

#include "scalar.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace regionTypes
{
    const dimensionedScalar kappa0_("kappa0", dimensionSet(-1, -3, 3, 0, 0, 2, 0), 116);
    const dimensionedScalar EKappa_("EKappa", dimensionSet(1, 2, -2, 0, -1, 0, 0), 15000);
    const dimensionedScalar ELambda_("ELambda", dimensionSet(1, 2, -2, 0, -1, 0, 0), 20000);
    const dimensionedScalar DLambda0_("DLambda0", dimensionSet(0, 2, -1, 0, 0, 0, 0), 1e-10);
    const dimensionedScalar fComp_("fComp", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.06);
    const dimensionedScalar HEC_("HEC", dimensionSet(1, 2, -2, 0, -1, 0, 0), 42000);
    const dimensionedScalar HAD_("HEC", dimensionSet(1, 2, -2, 0, -1, 0, 0), 42000);
    const dimensionedScalar aA_("aA", dimensionSet(0, 1, -1, 0, 0, 0, 0), 3.53e-5);
    const dimensionedScalar aD_("aD", dimensionSet(0, 1, -1, 0, 0, 0, 0), 1.42e-4);
    const dimensionedScalar TRefMu1_("TRefMu1", dimensionSet(0, 0, 0, 1, 0, 0, 0), 542.05);
    const dimensionedScalar TRefMu2_("TRefMu2", dimensionSet(0, 0, 0, 1, 0, 0, 0), 144.15);
    const dimensionedScalar TRefP1_("TRefP1", dimensionSet(0, 0, 0, 1, 0, 0, 0), 3816.44);
    const dimensionedScalar TRefP2_("TRefP2", dimensionSet(0, 0, 0, 1, 0, 0, 0), 46.13);
    const dimensionedScalar pDim_("pDim", dimensionSet(1, -1, -2, 0, 0, 0, 0), 1);
    const dimensionedScalar muDim_("muDim", dimensionSet(1,-1, -1, 0, 0, 0, 0), 1);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
