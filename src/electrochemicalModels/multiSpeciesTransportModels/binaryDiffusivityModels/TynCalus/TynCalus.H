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
    Foam::binaryLiquidDiffusivityModels::TynCalus

Description
   TynCalus binary diffusivity model.
    
SourceFiles
    TynCalus.C

\*---------------------------------------------------------------------------*/

#ifndef TynCalus_H
#define TynCalus_H

#include "binaryDiffusivityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// namespace binaryLiquidDiffusivityModels
// {

/*---------------------------------------------------------------------------*\
                           Class TynCalus Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class TynCalus
:
     public binaryDiffusivityModel<ThermoType>
{
    // Private data
          const dictionary binaryDiffDict_;
    
        //- Dynamic viscosity of the solvent
        dimensionedScalar mu2_;  
          
        // specie Molar volume at its normal boiling point
        scalar V12_;

        // Units correction factor
        dimensionedScalar unitsCorr_; 


public:

    //- Runtime type information
    TypeName("TynCalus");


    // Constructors

        //- construct from components
        TynCalus
        (
            const ThermoType& thermo,
            const word& name1,
            const word& name2,
            const IOdictionary& dict,
            const phaseModel& phase
        );


    // Destructor

        ~TynCalus()
        {}


         //- Binary diffusion coefficient [m^2/s]
         virtual tmp<volScalarField> D() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// } // End namespace binaryLiquidDiffusivityModels
} // End namespace Foam

#ifdef NoRepository
    #include "TynCalus.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
