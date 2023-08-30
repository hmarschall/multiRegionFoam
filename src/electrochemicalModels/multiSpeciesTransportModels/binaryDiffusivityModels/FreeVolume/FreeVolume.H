/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::binaryLiquidDiffusivityModels::FreeVolume

Description
   FreeVolume binary diffusivity model.
    
SourceFiles
    FreeVolume.C

\*---------------------------------------------------------------------------*/

#ifndef FreeVolume_H
#define FreeVolume_H

#include "binaryDiffusivityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// namespace binaryLiquidDiffusivityModels
// {

/*---------------------------------------------------------------------------*\
                           Class FreeVolume Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class FreeVolume
:
    public binaryDiffusivityModel<ThermoType>
{
    // Private data
    
        const dictionary binaryDiffDict_;
            
        // specie free volume
         //scalar vf1_;
        // specie denisty 
        dimensionedScalar rho2_;
       //- Boltzmann constant value (J/K)
        const dimensionedScalar Kb_;
       //- Specie diametre 1 
        //const scalar d1_;
        //- Specie diametre 2
         scalar d2_;
      // Gamma constant value
        const scalar Gamma_;
      //- Constant of the model
        scalar A_, B_; 
      // Units correction factor
        dimensionedScalar unitsCorr_;   

public:

    //- Runtime type information
    TypeName("FreeVolume");


    // Constructors

        //- construct from components
        FreeVolume
        (
            const ThermoType& thermo,
            const word& name1,
            const word& name2,
            const IOdictionary& dic,
            const phaseModel& phase       
        );


    // Destructor

        ~FreeVolume()
        {}


    // Member Functions

         //- Binary diffusion coefficient [m^2/s]
         virtual tmp<volScalarField> D() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// } // End namespace binaryLiquidDiffusivityModels
} // End namespace Foam

#ifdef NoRepository
    #include "FreeVolume.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
#endif

// ************************************************************************* //
