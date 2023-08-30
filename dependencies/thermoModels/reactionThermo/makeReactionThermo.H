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

InClass
    Foam::rhoReactionThermo

Description

\*---------------------------------------------------------------------------*/

// Does not work somehow
// typedef CThermo                                                    \
// <                                                                          \
//     BaseReactionThermo,                                                    \                                                                \
//     Mixture                                                                \
//     <                                                                      \
//         ThermoPhys                                                         \
//     >                                                                      \
// > CThermo##Mixture##ThermoPhys;                                            \
//                                                                            \
// defineTemplateTypeNameAndDebugWithName                                     \
// (                                                                          \
//     CThermo##Mixture##ThermoPhys,                                          \
//     (#CThermo"<"#Mixture"<"#ThermoPhys">>").c_str(),                       \
//     0                                                                      \
// );                                                                                                          

#ifndef makeReactionThermo_H
#define makeReactionThermo_H

#include "addToRunTimeSelectionTable.H"
#include "makeThermo.H"
#include "SpecieMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// defineTemplateTypeNameAndDebugWithName  definiert in Zeile 1 das ganze thermo object in "code"
//                                          definiert in Zeile 2 den Eintrag in "String" wonach der compiler suchen soll
//                                          d.h. hier
//  heRhoThermo<multiComponentMixture<sutherlandTransport<specieThermo<janafThermo<perfectGas<specie>>> 
// NOT heRhoThermo<rhoReactionThermo,multiComponentMixture<sutherlandTransport<specieThermo<janafThermo<perfectGas<specie>>>                         
// Bem: Ohne den rhoReactionThermo-Eintrag vor dem Komma
// Bsp. Of19
// hePsiThermo<homogeneousMixture<sutherland<janaf<perfectGas<specie>>,sensibleEnthalpy>>>

// this would do it as well     #CThermo"<"#Mixture"<"#ThermoPhys">>"
// but then the user would tell in thermoPhysicalProperties                                 \
// heRhoThermo<multiComponentMixture<gasThermoPhysics>>
// but one does not know for instance what the typedef gasThermoPhysics means

// #define defineThermoPhysicsReactionThermo(BaseReactionThermo,CThermo,Mixture,ThermoPhys) \
//                                                                                \
// typedef CThermo<BaseReactionThermo,Mixture<ThermoPhys>> CThermo##Mixture##ThermoPhys; \
//                                                                                   \
// defineTemplateTypeNameAndDebugWithName                                        \
// (                                                                             \
//     CThermo##Mixture##ThermoPhys,                                             \
//     (#CThermo"<" + Mixture<ThermoPhys>::typeName() + ">").c_str(),                                 \
//     0                                                                         \
// );

#define defineThermoPhysicsReactionThermo(BaseReactionThermo,CThermo,Mixture,ThermoPhys) \
                                                                               \
    typedef CThermo                                                            \
    <                                                                          \
        BaseReactionThermo,                                                    \
        SpecieMixture                                                          \
        <                                                                      \
            Mixture                                                            \
            <                                                                  \
                ThermoPhys                                                     \
            >                                                                  \
        >                                                                      \
    > CThermo##Mixture##ThermoPhys;                                            \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        CThermo##Mixture##ThermoPhys,                                          \
        (#CThermo"<" + Mixture<ThermoPhys>::typeName() + ">").c_str(),         \
        0                                                                      \
    )
    

#define makeThermoPhysicsReactionThermos(BaseThermo,BaseReactionThermo,CThermo,Mixture,ThermoPhys) \
                                                                               \
    defineThermoPhysicsReactionThermo                                          \
    (                                                                          \
        BaseReactionThermo,                                                    \
        CThermo,                                                               \
        Mixture,                                                               \
        ThermoPhys                                                             \
    );                                                                         \
                                                                               \
    addThermoPhysicsThermo(basicThermo, CThermo##Mixture##ThermoPhys);         \
    addThermoPhysicsThermo(fluidThermo, CThermo##Mixture##ThermoPhys);         \
    addThermoPhysicsThermo(BaseThermo, CThermo##Mixture##ThermoPhys);          \
    addThermoPhysicsThermo(BaseReactionThermo, CThermo##Mixture##ThermoPhys)


#define makeThermoPhysicsReactionThermo(BaseReactionThermo,CThermo,Mixture,ThermoPhys) \
                                                                               \
    defineThermoPhysicsReactionThermo                                          \
    (                                                                          \
        BaseReactionThermo,                                                    \
        CThermo,                                                               \
        Mixture,                                                               \
        ThermoPhys                                                             \
    );                                                                         \
                                                                               \
                                                                               \
    addThermoPhysicsThermo(BaseReactionThermo, CThermo##Mixture##ThermoPhys)


#define makeReactionThermos(BaseThermo,BaseReactionThermo,CThermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie) \
                                                                               \
    typedefThermoPhysics(Transport,Type,Thermo,EqnOfState,Specie);             \
                                                                               \
    makeThermoPhysicsReactionThermos                                           \
    (                                                                          \
        BaseThermo,                                                            \
        BaseReactionThermo,                                                    \
        CThermo,                                                               \
        Mixture,                                                               \
        Transport##Type##Thermo##EqnOfState##Specie                            \
    )

#define makeReactionThermo(BaseReactionThermo,CThermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie) \
                                                                               \
    typedefThermoPhysics(Transport,Type,Thermo,EqnOfState,Specie);             \
                                                                               \
    makeThermoPhysicsReactionThermo                                            \
    (                                                                          \
        BaseReactionThermo,                                                    \
        CThermo,                                                               \
        Mixture,                                                               \
        Transport##Type##Thermo##EqnOfState##Specie                            \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
