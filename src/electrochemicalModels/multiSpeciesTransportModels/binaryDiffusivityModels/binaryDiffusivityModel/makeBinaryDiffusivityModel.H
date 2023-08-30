/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

InClass
    Foam::basicPsiThermo

Description
    Macros for creating 'basic' compresibility-based thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeBinaryDiffusivityModel_H
#define makeBinaryDiffusivityModel_H

#include "addToRunTimeSelectionTable.H"
#include "binaryDiffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeBinaryDiffusivityModelTypes(binaryDiffusivityModel,Thermo)         \
                                                                            \
        typedef binaryDiffusivityModel<Thermo>                             \
            binaryDiffusivityModel##Thermo;                                \
                                                                                \
defineTemplateRunTimeSelectionTable(binaryDiffusivityModel##Thermo, dictionary); 


#define makeBinaryDiffusivityModel(BinaryDiffusivityModel,Thermo)                 \
                                                                              \
typedef BinaryDiffusivityModel<Thermo>                                                \
    BinaryDiffusivityModel##Thermo;                                                   \
                                                                              \
                                                                           \
                                                                            \
defineNamedTemplateTypeNameAndDebug                                        \
(                                                                             \
    BinaryDiffusivityModel##Thermo,                                                   \
    0                                                                         \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    binaryDiffusivityModel##Thermo,                                               \
    BinaryDiffusivityModel##Thermo,                                                   \
    dictionary                                                                    \
)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
