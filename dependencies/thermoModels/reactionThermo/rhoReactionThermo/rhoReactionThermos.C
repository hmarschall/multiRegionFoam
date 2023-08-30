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

#include "rhoReactionThermo.H"
#include "heRhoThermo.H"

#include "specie.H"
#include "perfectGas.H"
#include "incompressiblePerfectGas.H"
#include "hConstThermo.H"
#include "janafThermo.H"
#include "sensibleEnthalpy.H"
#include "specieThermo.H"
#include "rhoConst.H"
#include "rPolynomial.H"
#include "perfectFluid.H"
#include "adiabaticPerfectFluid.H"
#include "Boussinesq.H"

#include "constTransport.H"
#include "sutherlandTransport.H"
#include "WLFTransport.H"

#include "thermoPhysicsTypes.H"

#include "homogeneousMixture.H"
// #include "inhomogeneousMixture.H"
// #include "veryInhomogeneousMixture.H"
// // #include "dieselMixture.H"
#include "multiComponentMixture.H"
#include "reactingMixture.H"
#include "singleComponentMixture.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// makeReactionThermos
// (
//     // calls macro in makerhoReactionThermo
//     // which adds the objects to the runTime table
//     rhoThermo,              // BaseThermo
//     rhoReactionThermo,      // BaseReactionThermo 
//     heRhoThermo,            // CThermo
//     homogeneousMixture,     // Mixture
//     constTransport,         // Transport
//     sensibleEnthalpy,       // Type
//     hConstThermo,           // Thermo
//     perfectGas,             // EqnOfState
//     specie                  // Specie
// );

// makeReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     inhomogeneousMixture,
//     constTransport,
//     sensibleEnthalpy,
//     hConstThermo,
//     perfectGas,
//     specie
// );

// makeReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     veryInhomogeneousMixture,
//     constTransport,
//     sensibleEnthalpy,
//     hConstThermo,
//     perfectGas,
//     specie
// );

// makeReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     homogeneousMixture,
//     sutherlandTransport,
//     sensibleEnthalpy,
//     janafThermo,
//     perfectGas,
//     specie
// );

// makeReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     inhomogeneousMixture,
//     sutherlandTransport,
//     sensibleEnthalpy,
//     janafThermo,
//     perfectGas,
//     specie
// );

// makeReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     veryInhomogeneousMixture,
//     sutherlandTransport,
//     sensibleEnthalpy,
//     janafThermo,
//     perfectGas,
//     specie
// );


// makeReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     homogeneousMixture,
//     constTransport,
//     sensibleEnthalpy,
//     hConstThermo,
//     incompressiblePerfectGas,
//     specie
// );

// makeReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     inhomogeneousMixture,
//     constTransport,
//     sensibleEnthalpy,
//     hConstThermo,
//     incompressiblePerfectGas,
//     specie
// );

// makeReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     veryInhomogeneousMixture,
//     constTransport,
//     sensibleEnthalpy,
//     hConstThermo,
//     incompressiblePerfectGas,
//     specie
// );

// makeReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     homogeneousMixture,
//     sutherlandTransport,
//     sensibleEnthalpy,
//     janafThermo,
//     incompressiblePerfectGas,
//     specie
// );

// makeReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     inhomogeneousMixture,
//     sutherlandTransport,
//     sensibleEnthalpy,
//     janafThermo,
//     incompressiblePerfectGas,
//     specie
// );

// makeReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     veryInhomogeneousMixture,
//     sutherlandTransport,
//     sensibleEnthalpy,
//     janafThermo,
//     incompressiblePerfectGas,
//     specie
// );

// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Single-component thermo for internal energy

makeThermoPhysicsReactionThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    constGasEThermoPhysics
);

makeThermoPhysicsReactionThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    gasEThermoPhysics
);

makeThermoPhysicsReactionThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    constIncompressibleGasEThermoPhysics
);

makeThermoPhysicsReactionThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    incompressibleGasEThermoPhysics
);

// makeThermoPhysicsReactionThermo
// (
//     rhoReactionThermo,
//     heRhoThermo,
//     singleComponentMixture,
//     icoPoly8EThermoPhysics
// );

makeThermoPhysicsReactionThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    constFluidEThermoPhysics
);

makeThermoPhysicsReactionThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    constAdiabaticFluidEThermoPhysics
);

makeThermoPhysicsReactionThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    constEThermoPhysics
);

makeReactionThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
    Boussinesq,
    specie
);

makeReactionThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    janafThermo,
    Boussinesq,
    specie
);

makeReactionThermo 
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    WLFTransport,
    sensibleInternalEnergy,
    eConstThermo,
    rhoConst,
    specie
);

// Multi-component thermo for sensible enthalpy

makeThermoPhysicsReactionThermos        
(
    rhoThermo,              // BaseThermo
    rhoReactionThermo,      // BaseReactionThermo
    heRhoThermo,            // CThermo
    multiComponentMixture,  // Mixture
    constGasHThermoPhysics  // ThermoPhys
);

makeThermoPhysicsReactionThermos        
(
    rhoThermo,              // BaseThermo
    rhoReactionThermo,      // BaseReactionThermo
    heRhoThermo,            // CThermo
    multiComponentMixture,  // Mixture
    gasEThermoPhysics  // ThermoPhys
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    gasHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constIncompressibleGasHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    incompressibleGasHThermoPhysics
);

// makeThermoPhysicsReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     multiComponentMixture,
//     icoPoly8HThermoPhysics
// );

// makeThermoPhysicsReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     multiComponentMixture,
//     icoPoly8TranspJanafHThermoPhysics
// );

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constFluidHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constAdiabaticFluidHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constHThermoPhysics
);

// // Reaction thermo for internal energy

// makeThermoPhysicsReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     reactingMixture,
//     constGasHThermoPhysics
// );

// makeThermoPhysicsReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     reactingMixture,
//     gasEThermoPhysics
// );

// makeThermoPhysicsReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     reactingMixture,
//     constIncompressibleGasEThermoPhysics
// );

// makeThermoPhysicsReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     reactingMixture,
//     incompressibleGasEThermoPhysics
// );

// // makeThermoPhysicsReactionThermos
// // (
// //     rhoThermo,
// //     rhoReactionThermo,
// //     heRhoThermo,
// //     reactingMixture,
// //     icoPoly8EThermoPhysics
// // );

// makeThermoPhysicsReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     reactingMixture,
//     constFluidEThermoPhysics
// );

// makeThermoPhysicsReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     reactingMixture,
//     constAdiabaticFluidEThermoPhysics
// );

// makeThermoPhysicsReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     reactingMixture,
//     constEThermoPhysics
// );




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
