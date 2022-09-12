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

#include "capillaryInterface.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionInterfaces
{
    defineTypeNameAndDebug(capillaryInterface, 0);

    addToRunTimeSelectionTable
    (
        regionInterface, 
        capillaryInterface, 
        IOdictionary 
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
    
Foam::regionInterfaces::capillaryInterface::capillaryInterface
( 
    const Time& runTime,   
    const fvPatch& patchA, 
    const fvPatch& patchB  
)
:
    regionInterface(runTime, patchA, patchB),

    sigma0_
    (
        interfaceProperties().subDict(name()).lookup("sigma")
    ),
    sigma_
    (
        areaScalarField
        (
            IOobject
            (
                "sigma" + aMesh().mesh().name() + patchA.name(),
                runTime.timeName(),
                aMesh().thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh(),
            sigma0_,
            zeroGradientFaPatchScalarField::typeName
        )
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionInterfaces::capillaryInterface::correct()
{
    //TO-Do call function to calculate new sigma_ with interface equation of state 
    // for contaminated surfaces 
}

Foam::scalar Foam::regionInterfaces::capillaryInterface::getMinDeltaT()
{
    scalar minDeltaT = GREAT;

    if (meshA().foundObject<volScalarField>("rho"))
    {
        const scalarField& dE = aMesh().lPN();
        scalar minDE = gMin(dE);
        const scalarField& rhoA = meshA().lookupObject<volScalarField>("rho");
        scalar minRhoA = gMin(rhoA);

        minDeltaT =
            0.9*
            sqrt
            (
                minRhoA*minDE*minDE*minDE/
                2.0/M_PI/(sigma0_.value() + SMALL)
            );

    }

    return minDeltaT;
}


// ************************************************************************* //
