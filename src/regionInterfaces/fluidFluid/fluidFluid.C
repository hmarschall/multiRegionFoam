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

#include "fluidFluid.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionInterfaces
{
    defineTypeNameAndDebug(fluidFluid, 0);

    addToRunTimeSelectionTable
    (
        regionInterface, 
        fluidFluid, 
        IOdictionary 
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
    
Foam::regionInterfaces::fluidFluid::fluidFluid
( 
    const Time& runTime,   
    const fvPatch& patchA, 
    const fvPatch& patchB  
)
:
    regionInterface(runTime, patchA, patchB),
    
    runTime_(runTime),
    
    patchA_(patchA),
    patchB_(patchB),
    
    meshA_(patchA.boundaryMesh().mesh()), 
    meshB_(patchB.boundaryMesh().mesh()), 
    
    U_
    (
        meshA().lookupObject<volVectorField>("U") 
    ),    
    phi_
    (
        meshA().lookupObject<surfaceScalarField>("phi")
    ),    
    rhoA_
    (
        this->subDict(meshA().name()).lookup("rho")
    ),
    rhoB_
    (
        this->subDict(meshB().name()).lookup("rho")
    ),
    muA_
    (
        this->subDict(meshA().name()).lookup("mu")
    ),
    muB_
    (
        this->subDict(meshB().name()).lookup("mu")
    ),
    sigma0_
    (
        this->lookup("cleanSurfaceTension")
    ),
    g_
    (
        this->lookup("g")
    )
       
{
    // add
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionInterfaces::fluidFluid::~fluidFluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// virtual functions from regionInterface.H 
void Foam::regionInterfaces::fluidFluid::attach() 
{
    // do nothing, add as required
}

void Foam::regionInterfaces::fluidFluid::detach() 
{
    // do nothing, add as required
}

void Foam::regionInterfaces::fluidFluid::updateInterpolatorAndGlobalPatches() 
{
    // do nothing, add as required
}

// ************************************************************************* //
