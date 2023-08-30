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

Typedefs
    Foam::thermoPhysicsTypes

Description
    Type definitions for thermo-physics models

\*---------------------------------------------------------------------------*/

#ifndef thermoPhysicsTypes_H
#define thermoPhysicsTypes_H

#include "specie.H"

// Equation of state
#include "perfectGas.H"
#include "icoPolynomial.H"
#include "perfectFluid.H"
#include "rhoConst.H"
#include "adiabaticPerfectFluid.H"
#include "incompressiblePerfectGas.H"
#include "rPolynomial.H"
#include "PengRobinsonGas.H"
#include "Boussinesq.H"

// Thermos
#include "hConstThermo.H"
#include "eConstThermo.H"
#include "hPolynomialThermo.H"
#include "janafThermo.H"
#include "specieThermo.H"
#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"

// Transport


#include "sutherlandTransport.H"
#include "constTransport.H"
#include "polynomialTransport.H"




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // typedef
    // sutherlandTransport
    // <
    //     species::specieThermo
    //     <
    //         hConstThermo // janafThermo // hConstThermo
    //         <
    //             // It might be a good idea to introduce the specie into the template
    //             // >>> Since then one can change between the
    //             // possiblity to create a charged specie (chargedSpecie) and normal specie or having only simple species
    //             perfectGas<specie>   //perfectGas   // perfectGas<specie>      
    //         >,
    //         sensibleEnthalpy
    //     >
    // > gasThermoPhysics;

    // thermo physics types based on sensibleEnthalpy
    typedef
    constTransport
    <
        species::specieThermo
        <
            hConstThermo
            <
                perfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > constGasHThermoPhysics;

    typedef
    sutherlandTransport
    <
        species::specieThermo
        <
            janafThermo
            <
                perfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > gasHThermoPhysics;

    typedef
    constTransport
    <
        species::specieThermo
        <
            hConstThermo
            <
                incompressiblePerfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > constIncompressibleGasHThermoPhysics;

    typedef
    sutherlandTransport
    <
        species::specieThermo
        <
            janafThermo
            <
                incompressiblePerfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > incompressibleGasHThermoPhysics;

    typedef
    polynomialTransport
    <
        species::specieThermo
        <
            hPolynomialThermo
            <
                icoPolynomial<specie, 8>,
                8
            >,
            sensibleEnthalpy
        >,
        8
    > icoPoly8HThermoPhysics;

    typedef
    polynomialTransport
    <
        species::specieThermo
        <
            janafThermo
            <
                perfectGas<specie>
            >,
            sensibleEnthalpy
        >,
        8
    > icoPoly8TranspJanafHThermoPhysics;

    typedef
    constTransport
    <
        species::specieThermo
        <
            hConstThermo
            <
                perfectFluid<specie>
            >,
            sensibleEnthalpy
        >
    >
    constFluidHThermoPhysics;

    typedef
    constTransport
    <
        species::specieThermo
        <
            hConstThermo
            <
                rPolynomial<specie>
            >,
            sensibleEnthalpy
        >
    >
    constrPolFluidHThermoPhysics;

    typedef
    constTransport
    <
        species::specieThermo
        <
            hConstThermo
            <
                adiabaticPerfectFluid<specie>
            >,
            sensibleEnthalpy
        >
    >
    constAdiabaticFluidHThermoPhysics;

    typedef
    constTransport
    <
        species::specieThermo
        <
            hConstThermo
            <
                rhoConst<specie>
            >,
            sensibleEnthalpy
        >
    >
    constHThermoPhysics;


    // thermo physics types based on sensibleInternalEnergy
    typedef
    constTransport
    <
        species::specieThermo
        <
            eConstThermo
            <
                perfectGas<specie>
            >,
            sensibleInternalEnergy
        >
    > constGasEThermoPhysics;

    typedef
    sutherlandTransport
    <
        species::specieThermo
        <
            janafThermo
            <
                perfectGas<specie>
            >,
            sensibleInternalEnergy
        >
    > gasEThermoPhysics;

    typedef
    constTransport
    <
        species::specieThermo
        <
            eConstThermo
            <
                incompressiblePerfectGas<specie>
            >,
            sensibleInternalEnergy
        >
    > constIncompressibleGasEThermoPhysics;

    typedef
    sutherlandTransport
    <
        species::specieThermo
        <
            janafThermo
            <
                incompressiblePerfectGas<specie>
            >,
            sensibleInternalEnergy
        >
    > incompressibleGasEThermoPhysics;

    typedef
    polynomialTransport
    <
        species::specieThermo
        <
            hPolynomialThermo
            <
                icoPolynomial<specie, 8>,
                8
            >,
            sensibleInternalEnergy
        >,
        8
    > icoPoly8EThermoPhysics;

    typedef
    polynomialTransport
    <
        species::specieThermo
        <
            janafThermo
            <
                perfectGas<specie>
            >,
            sensibleInternalEnergy
        >,
        8
    > icoPoly8TranspJanafEThermoPhysics;

    typedef
    constTransport
    <
        species::specieThermo
        <
            eConstThermo
            <
                perfectFluid<specie>
            >,
            sensibleInternalEnergy
        >
    >
    constFluidEThermoPhysics;

    typedef
    constTransport
    <
        species::specieThermo
        <
            eConstThermo
            <
                perfectFluid<specie>
            >,
            sensibleInternalEnergy
        >
    >
    constrPolFluidEThermoPhysics;

    typedef
    constTransport
    <
        species::specieThermo
        <
            eConstThermo
            <
                adiabaticPerfectFluid<specie>
            >,
            sensibleInternalEnergy
        >
    >
    constAdiabaticFluidEThermoPhysics;

    typedef
    constTransport
    <
        species::specieThermo
        <
            eConstThermo
            <
                rhoConst<specie>
            >,
            sensibleInternalEnergy
        >
    >
    constEThermoPhysics;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

