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

#include "makeReactionThermo.H"
#include "thermoPhysicsTypes.H"
// #include "solidThermoPhysicsTypes.H"

#include "chemistryReader.H"
#include "foamChemistryReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// // Fluid chemistry readers based on sensibleEnthalpy

// makeChemistryReader(constGasHThermoPhysics);
// // makeChemistryReader(gasHThermoPhysics);
// // makeChemistryReader(constIncompressibleGasHThermoPhysics);
// // makeChemistryReader(incompressibleGasHThermoPhysics);
// // makeChemistryReader(icoPoly8HThermoPhysics);
// // makeChemistryReader(constFluidHThermoPhysics);
// // makeChemistryReader(constAdiabaticFluidHThermoPhysics);
// // makeChemistryReader(constHThermoPhysics);


// makeChemistryReaderType(foamChemistryReader, constGasHThermoPhysics);
// // makeChemistryReaderType(foamChemistryReader, gasHThermoPhysics);
// // makeChemistryReaderType
// // (
// //     foamChemistryReader,
// //     constIncompressibleGasHThermoPhysics
// // );
// // makeChemistryReaderType(foamChemistryReader, incompressibleGasHThermoPhysics);
// // makeChemistryReaderType(foamChemistryReader, icoPoly8HThermoPhysics);
// // makeChemistryReaderType(foamChemistryReader, constFluidHThermoPhysics);
// // makeChemistryReaderType(foamChemistryReader, constAdiabaticFluidHThermoPhysics);
// // makeChemistryReaderType(foamChemistryReader, constHThermoPhysics);



// // // Fluid chemistry readers based on sensibleInternalEnergy

// // makeChemistryReader(constGasEThermoPhysics);
// // makeChemistryReader(gasEThermoPhysics);
// // makeChemistryReader(constIncompressibleGasEThermoPhysics);
// // makeChemistryReader(incompressibleGasEThermoPhysics);
// // makeChemistryReader(icoPoly8EThermoPhysics);
// // makeChemistryReader(constFluidEThermoPhysics);
// // makeChemistryReader(constAdiabaticFluidEThermoPhysics);
// // makeChemistryReader(constEThermoPhysics);


// // makeChemistryReaderType(foamChemistryReader, constGasEThermoPhysics);
// // makeChemistryReaderType(foamChemistryReader, gasEThermoPhysics);
// // makeChemistryReaderType
// // (
// //     foamChemistryReader,
// //     constIncompressibleGasEThermoPhysics
// // );
// // makeChemistryReaderType(foamChemistryReader, incompressibleGasEThermoPhysics);
// // makeChemistryReaderType(foamChemistryReader, icoPoly8EThermoPhysics);
// // makeChemistryReaderType(foamChemistryReader, constFluidEThermoPhysics);
// // makeChemistryReaderType(foamChemistryReader, constAdiabaticFluidEThermoPhysics);
// // makeChemistryReaderType(foamChemistryReader, constEThermoPhysics);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
