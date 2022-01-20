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

    sigma_
    (
        interfaceProperties().subDict(name()).lookup("sigma")
    ),
    sigmaPtr_(NULL)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<areaScalarField>
Foam::regionInterfaces::capillaryInterface::sigma() const
{         
    if
    (
        //meshA().foundObject<areaScalarField>("sigma")
        runTime().foundObject<areaScalarField>("sigma")
    )
    {
        // contaminated interface
        
        //sigmaPtr_ = const_cast<areaScalarField*>
        //    (&meshA().lookupObject<areaScalarField>("sigma"));
            
        sigmaPtr_.set
        (
            new areaScalarField
            (
               //meshA().lookupObject<areaScalarField>("sigma")
               runTime().lookupObject<areaScalarField>("sigma")
            )
        );   
    }
    else
    {
       // clean interface
       
       // sigmaPtr_ = new areaScalarField
       // (
       //     IOobject
       //     (
       //         "sigma",
       //         this->db().time().timeName(),
       //         aMesh().thisDb(),
       //         IOobject::NO_READ,
       //         IOobject::NO_WRITE
       //     ),
       //     aMesh(),
       //     sigma_,
       //     zeroGradientFaPatchVectorField::typeName
       // );
        
        
        sigmaPtr_.set
        (
            new areaScalarField
            (
                IOobject
                (
                    "sigma",
                    runTime().timeName(),
                    aMesh().thisDb(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                aMesh(),
                sigma_,
                zeroGradientFaPatchVectorField::typeName
            )
        );
    }      
    
    //return sigmaPtr;
} 


// ************************************************************************* //
