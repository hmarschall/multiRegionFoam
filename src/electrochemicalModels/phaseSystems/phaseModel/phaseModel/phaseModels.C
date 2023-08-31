/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"

#include "rhoThermo.H"
#include "rhoReactionThermo.H"

#include "phaseModel.H"
#include "ThermoPhaseModel.H"
#include "AnisothermalPhaseModel.H"
#include "MultiComponentPhaseModel.H"
#include "SurfaceMultiComponentPhaseModel.H"
#include "SurfaceElectrochemicalReactingPhaseModel.H"
#include "InertPhaseModel.H"
#include "MovingPhaseModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Volumetric electrochemistry
    // typedef
    //     // AnisothermalPhaseModel
    //     // <
    //         MultiComponentFickPhaseModel
    //         <
    //             ElectrochemicalReactingPhaseModel
    //             <
    //                 MovingPhaseModel
    //                 <
    //                     ThermoPhaseModel<phaseModel, rhoReactionThermo>
    //                 >
    //             >
    //         >
    //     // >
    //     reactingFickPhaseModel;

    // addNamedToRunTimeSelectionTable
    // (
    //     // (baseType,thisType,argNames,lookup)
    //     phaseModel,
    //     reactingFickPhaseModel,
    //     // phaseSystem,     // phaseSystem was/is an IOdictionary
    //     dictionary,
    //     reactingFickPhaseModel
    // );



    typedef
        AnisothermalPhaseModel
        <
            MultiComponentPhaseModel
            <
                InertPhaseModel
                <
                    MovingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoReactionThermo>
                    >
                >
            >
        >
        multiComponentPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        multiComponentPhaseModel,
        dictionary,
        multiComponentPhaseModel
    );


    // Surface coupled/ electrochemistry phase models
    // Reaction
    typedef
        AnisothermalPhaseModel
        <
            SurfaceMultiComponentPhaseModel
            <
                SurfaceElectrochemicalReactingPhaseModel
                <
                    MovingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoReactionThermo>
                    >
                >
            >
        >
        surfaceReactingFickPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        surfaceReactingFickPhaseModel,
        dictionary,
        surfaceReactingFickPhaseModel
    );

}

// ************************************************************************* //
