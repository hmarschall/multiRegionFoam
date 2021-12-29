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

    transportPropertiesA_
    (
        IOobject
        (
            "transportProperties",
            fileName(runTime.caseConstant()/meshA().name()),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    transportPropertiesB_
    (
        IOobject
        (
            "transportProperties",
            fileName(runTime.caseConstant()/meshB().name()),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    gravitationalProperties_
    (
        IOobject
        (
            "g",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

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
        transportPropertiesA_.lookup("rho")
    ),
    rhoB_
    (
        transportPropertiesB_.lookup("rho")
    ),
    muA_
    (
        transportPropertiesA_.lookup("mu")
    ),
    muB_
    (
        transportPropertiesB_.lookup("mu")
    ),
    sigma0_
    (
        interfaceProperties().subDict(name()).lookup("sigma")
    ),
    g_
    (
        gravitationalProperties_.lookup("g")
    )
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //





// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //





// ************************************************************************* //
