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
    Foam::Fickian

Description
    Fickian multispecies mass transport model.
  
SourceFiles
    Fickian.C

\*---------------------------------------------------------------------------*/

#ifndef Fickian_H
#define Fickian_H

#include "MultiSpeciesTransportModel.H"
#include "fvm.H"

#include "phaseSystem.H"
#include "phaseModel.H"
#include "basicSpecieMixture.H"
#include "binaryDiffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace multiSpeciesTransportModels
{
  
/*---------------------------------------------------------------------------*\
                             Class Fickian Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class Fickian
:
    public MultiSpeciesTransportModel<Thermo>
{  

private:
  
    // Protected data
    
        //- Access to molar fractions
        const PtrList<volScalarField>& X_;

        //- Multispecies diffusion coefficients
        PtrList<volScalarField> Dm_;

        //- Binary diffusivity models 
        PtrList<binaryDiffusivityModel<Thermo>> DijModels_; 

        //- Binary diffusion coefficients
        PtrList<volScalarField> Dij_;  

public:

    //- Runtime type information
    TypeName("Fickian");

    
    // Constructors

        //- Construct from components
        Fickian
        (
            const phaseModel& phase
        );


      //- Destructor
      virtual ~Fickian();


    // Member Functions

    //- Effective mass diffusivity
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const;

        //- Effective mass diffusivity
        //  for a given specie mass-fraction for patch [kg/m/s]
        virtual tmp<scalarField> DEff
        (
            const volScalarField& Yi,
            const label patchi
        ) const;

        //- Return the heat flux [W/m^2]
        virtual tmp<surfaceScalarField> q() const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

        //- Return the specie flux for the given specie mass-fraction [kg/m^2/s]
        virtual tmp<surfaceScalarField> j(const volScalarField& Yi) const;

        //- Return the source term for the given specie mass-fraction equation
        virtual tmp<fvScalarMatrix> divj(volScalarField& Yi) const;

        //- Update the diffusion coefficients
        virtual void correct();

        //- Return the const binary diffusion coefficient field
        const volScalarField& Dij(label i, label j) const;

        void updateBinaryDiffusionCoeff();

  };

  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
}   // End namespace multiSpeciesTransportModels
}	// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "Fickian.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
