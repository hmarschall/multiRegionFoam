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

\*---------------------------------------------------------------------------*/

#include "makeBinaryDiffusivityModel.H"
#include "binaryDiffusivityModel.H"
#include "rhoReactionThermo.H"
#include "rhoThermo.H"

#include "Fuller.H"
// #include "Wilke.H"
// #include "ChapmanEnskog.H"


// #include "StokeEinstein.H"
// #include "WilkeChang.H"
// #include "Wilke.H"
// #include "constantBinaryDiff.H"
// #include "FreeVolume.H"
// #include "FvPlusSE.H"
// #include "HaydukLaudie.H"
// #include "TynCalus.H"
// #include "Knudsen.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

//TODO:
// For now binaryDiffModel is also a template
// Probably need to be changed
// -> Define it here for rhoReactionThermo
makeBinaryDiffusivityModelTypes(binaryDiffusivityModel, rhoReactionThermo);

//- Gas binaryDiffusionModels
makeBinaryDiffusivityModel(Fuller, rhoReactionThermo);
// makeBinaryDiffusivityModel(ChapmanEnskog, rhoReactionThermo);
// makeBinaryDiffusivityModel(Wilke, rhoReactionThermo);
// makeBinaryDiffusivityModel(constantBinaryDiff, rhoReactionThermo);
// makeBinaryDiffusivityModel(Knudsen, rhoReactionThermo);

// //- Liquid binaryDiffusionModels
// makeBinaryDiffusivityModel(StokeEinstein, rhoReactionThermo);
// makeBinaryDiffusivityModel(WilkeChang, rhoReactionThermo);
// makeBinaryDiffusivityModel(FreeVolume, rhoReactionThermo);
// makeBinaryDiffusivityModel(FvPlusSE, rhoReactionThermo);
// makeBinaryDiffusivityModel(TynCalus, rhoReactionThermo);
// makeBinaryDiffusivityModel(HaydukLaudie, rhoReactionThermo);



// // rhoThermo
// makeBinaryDiffusivityModelTypes(binaryDiffusivityModel, rhoThermo);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
