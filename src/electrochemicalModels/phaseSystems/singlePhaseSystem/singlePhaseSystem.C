/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
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
#include "singlePhaseSystem.H"
#include "addToRunTimeSelectionTable.H"

#include "fvCFD.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(singlePhaseSystem, 0);
    addToRunTimeSelectionTable(phaseSystem, singlePhaseSystem, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singlePhaseSystem::singlePhaseSystem
(
    const fvMesh& mesh
)
:
    phaseSystem(mesh),

    phase_(phaseModels_[0]),
    surfaceElectrochemistry_(this->lookupOrDefault<Switch>("surfaceElectrochemistry", false))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::singlePhaseSystem::~singlePhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Foam::tmp<Foam::surfaceScalarField>
// Foam::singlePhaseSystem::phirMag() const
// {
//     return mag(phase_.phi());
// }


// Foam::tmp<Foam::surfaceScalarField>
// Foam::singlePhaseSystem::phiMagMax() const
// {
//     return mag(phase_.phi());
// }


void Foam::singlePhaseSystem::solve()
{
    Info << "Solve for single phase flow:" << endl;

    fvMesh& mesh = const_cast<fvMesh&>(mesh_);

    const Time& runTime = mesh.time();

    #include "createFields.H"

    Switch faceMomentum
    (
        pimple_.dict().lookupOrDefault<Switch>("faceMomentum", false)
    );

    Switch Y
    (
        pimple_.dict().lookupOrDefault<Switch>("Y", true)
    );

    // --- Pressure-velocity PIMPLE corrector loop
    while (pimple_.loop())
    {
        correct();

        #include "YEqns.H"

        Info << "Test debug " << endl;

        // if (faceMomentum)
        // {
        //     #include "pUf/UEqn.H"
        //     #include "pUf/pEqn.H"
        // }
    // //     else
    // //     {
            #include "pU/UEqn.H"
            #include "pU/pEqn.H"

            // #include "EEqn.H"
            // #include "TEqn.H"
    // //     }

        // correctKinematics();

    // //     if (pimple_.turbCorr())
    // //     {
    // //         correctTurbulence();
    // //     }
    }

    // Update the boundary conditions for the surface coupled cases
    if(surfaceElectrochemistry_)
    {
        forAll(phaseModels_, i)
        {
            phaseModels_[i].correctBC();
        }
    }
}


void Foam::singlePhaseSystem::correct()
{
    phase_.correct();
}

Foam::tmp<Foam::fvScalarMatrix> Foam::singlePhaseSystem::TEqn()
{

    // const volScalarField& alpha = phase_;

    // rhoThermo& thermo = phase_.thermoRef();
    // volScalarField& rho = thermo.rho();

    // volVectorField& U = phase_.URef();

    // if (!phase_.isothermal())
    // {
        tmp<fvScalarMatrix> tTEqn
        (
            phase_.TEqn()
        //  ==
        //    alpha*rho*(U&g)
        );
    // } 

    return tTEqn;
}

// ************************************************************************* //
