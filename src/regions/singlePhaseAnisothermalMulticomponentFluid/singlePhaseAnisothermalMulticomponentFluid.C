/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "fvCFD.H"

#include "singlePhaseAnisothermalMulticomponentFluid.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(singlePhaseAnisothermalMulticomponentFluid, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        singlePhaseAnisothermalMulticomponentFluid,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::singlePhaseAnisothermalMulticomponentFluid::singlePhaseAnisothermalMulticomponentFluid
(
    const Time& runTime,
    const word& regionName
)
:
    regionType(runTime, regionName),
    regionName_(regionName)
{
    phases_ = phaseSystem::New
    (
        this->mesh()
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::singlePhaseAnisothermalMulticomponentFluid::~singlePhaseAnisothermalMulticomponentFluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::singlePhaseAnisothermalMulticomponentFluid::correct()
{
    phases_->correctEnergyTransport();
    phases_->correctThermo();
    phases_->correct();
}

Foam::scalar Foam::regionTypes::singlePhaseAnisothermalMulticomponentFluid::getMinDeltaT()
{
    return GREAT;
}



void Foam::regionTypes::singlePhaseAnisothermalMulticomponentFluid::setCoupledEqns()
{
    word matrixSystemNameRegion =
    (
        phases_->phases()[0].thermo().T().name()
      + phases_->mesh().name() + "Mesh"
      + singlePhaseAnisothermalMulticomponentFluid::typeName + "Type"
      + "Eqn"
    );

    TEqn = phases_->TEqn();

    fvScalarMatrices.set
    (
        phases_->phases()[0].thermo().T().name()
      + phases_->mesh().name() + "Mesh"
      + singlePhaseAnisothermalMulticomponentFluid::typeName + "Type"
      + "Eqn",
        &TEqn()
    );
}

void Foam::regionTypes::singlePhaseAnisothermalMulticomponentFluid::solveRegion()
{
    phases_->solve();
}

void Foam::regionTypes::singlePhaseAnisothermalMulticomponentFluid::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::singlePhaseAnisothermalMulticomponentFluid::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::singlePhaseAnisothermalMulticomponentFluid::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::singlePhaseAnisothermalMulticomponentFluid::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::singlePhaseAnisothermalMulticomponentFluid::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //
