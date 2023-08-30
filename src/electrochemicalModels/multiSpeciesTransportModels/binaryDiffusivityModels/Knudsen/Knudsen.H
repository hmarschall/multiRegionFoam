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
    Foam::Knudsen

Description
    Knudsen diffusivity model.
    
    Knudsen diffusivity coefficient is obtained by
    @f[
        D_{ik} = r \cdot \frac{2}{3} \cdot \sqrt{\frac{8 R T}{\pi M_i}}
    @f]

SourceFiles
    Knudsen.C

\*---------------------------------------------------------------------------*/

#ifndef Knudsen_H
#define Knudsen_H

#include "binaryDiffusivityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// namespace binaryDiffusivityModels
// {

/*---------------------------------------------------------------------------*\
                           Class Knudsen Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class Knudsen
:
    public binaryDiffusivityModel<ThermoType>
{    
    // Private data
    
        const dictionary binaryDiffDict_;

        //- Perfect gas constant value (J/kmol/K)
        const scalar RR_, PI_;

        //- Specie molar weight
        scalar W1_;

        //- Mean pore diametre
        const dimensionedScalar dp_;
        
         // Units correction factor
        dimensionedScalar unitsCorr_; 
    

public:

    //- Runtime type information
    TypeName("Knudsen");


    // Constructors

        //- construct from components
        Knudsen
        (
            const ThermoType& thermo,
            const word name1,
            const word name2,            
            const IOdictionary& dict,
            const phaseModel& phase
        );


    // Destructor

        ~Knudsen()
        {}


    // Member Functions
        
        //- Knudsen diffusion coefficient in m^2/s
        virtual tmp<volScalarField> D() const;                
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
    #include "Knudsen.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
