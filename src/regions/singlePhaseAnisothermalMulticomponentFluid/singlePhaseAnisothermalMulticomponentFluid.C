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
    const fvMesh& mesh,
    const word& regionName
)
:
    regionType(mesh, regionName),
    mesh_(mesh),  
    regionName_(regionName),
    dict_
    (
        IOobject
        (
            "phaseProperties",
            this->time().constant(),
            *this,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    electrochemicalPhaseModel_(phaseModel::New(dict_, *this)),
    pimple_(*this),
    p_rgh_
    (
        IOobject
        (
            "p_rgh",
            this->time().timeName(),
            *this,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        *this
    )
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::singlePhaseAnisothermalMulticomponentFluid::~singlePhaseAnisothermalMulticomponentFluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::singlePhaseAnisothermalMulticomponentFluid::correct()
{
    Info << "Correct electrochemistry phase model " << endl;
    electrochemicalPhaseModel_->correctEnergyTransport();
    electrochemicalPhaseModel_->correctThermo();
    electrochemicalPhaseModel_->correct();
}


void Foam::regionTypes::singlePhaseAnisothermalMulticomponentFluid::setRDeltaT()
{
    // do nothing, add as required
}

void Foam::regionTypes::singlePhaseAnisothermalMulticomponentFluid::updateFields()
{
    // do nothing, add as required
}


void Foam::regionTypes::singlePhaseAnisothermalMulticomponentFluid::setCoupledEqns()
{
    
// TODO: pahse_.setCoupleUEqn();
//       phase_.setEeqn();

}

void Foam::regionTypes::singlePhaseAnisothermalMulticomponentFluid::solveRegion()
{

    Info << "Solve for single phase flow:" << endl;

    fvMesh& mesh = const_cast<fvMesh&>(mesh_);

    const Time& runTime = mesh.time();

    // #include "createFields.H"
    #include "createFields.H"

    Switch faceMomentum
    (
        pimple_.dict().lookupOrDefault<Switch>("faceMomentum", false)
    );

        // --- Pressure-velocity PIMPLE corrector loop
    while (pimple_.loop())
    {
        electrochemicalPhaseModel_->correct();

        #include "YEqns.H"

        // if (faceMomentum)
        // {
        //     #include "pUf/UEqn.H"
        //     #include "pUf/pEqn.H"
        // }
        // else
        // {
            #include "pU/UEqn.H"
            #include "pU/pEqn.H"
        // }

        // correctKinematics();

        // if (pimple_.turbCorr())
        // {
        //     correctTurbulence();
        // }
    }

}

// ************************************************************************* //
