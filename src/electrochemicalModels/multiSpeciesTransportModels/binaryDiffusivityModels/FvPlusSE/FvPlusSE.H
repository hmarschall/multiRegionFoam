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
    Foam::FvPlusSE

Description
    FvPlusSE diffusivity model.
    

SourceFiles
    FvPlusSE.C

\*---------------------------------------------------------------------------*/

#ifndef FvPlusSE_H
#define FvPlusSE_H

#include "binaryDiffusivityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
//namespace binaryLiquidDiffusivityModels
//{
/*---------------------------------------------------------------------------*\
                           Class FvPlusStokeEi Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class FvPlusSE
:
    public binaryDiffusivityModel<ThermoType>
{
    // Private data

        const dictionary binaryDiffDict_;    

        // Boltzmann constant value (J/K)
        const scalar Kb_;

        //- Dynamic viscosity of the solvent
        const scalar mu_;

        //- Density of the solvent
        const scalar rho1_;   
        
        //- Density of the solute    
        const scalar rho2_;
        
        // Pi
        const scalar PI_;

        // Equation Constants
        scalar A_, B_, C_,Gamma_ ,b_;
        
        //Self Diffusion Coeffs from FreeVolume and StokeEinstien
        scalar Dself_Fv,Dself_SE;
        
        // Units correction factor
        dimensionedScalar unitsCorr_;     

public:

    //- Runtime type information
    TypeName("FvPlusSE");


    // Constructors

        //- construct from components
        FvPlusSE
        (
            const ThermoType& thermo,
            const word& name1,
            const word& name2,
            const IOdictionary& dict,
            const phaseModel& phase
        );


    // Destructor

        ~FvPlusSE()
        {}


    // Member Functions

        //- FvPlusSE diffusion coefficient in m^2/s
        virtual tmp<volScalarField> D() const;                
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// } // End namespace binaryLiquidDiffusivityModels
} // End namespace Foam
#ifdef NoRepository
    #include "FvPlusSE.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
