/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

Class
    Foam::turbulenceThermophysicalTransportModels::FickianBosanquet

Description
    Base class for multi-component FickianBosanquet based temperature gradient heat
    flux models with optional Soret thermal diffusion of species.

    The mixture diffusion coefficients are specified as Function2<scalar>s of
    pressure and temperature but independent of composition.

    The heat flux source is implemented as an implicit energy correction to the
    temperature gradient based flux source.  At convergence the energy
    correction is 0.

SourceFiles
    FickianBosanquet.C

\*---------------------------------------------------------------------------*/

// #include "Function2.H"
// TODO: Here it need access to the binaryDiffusion classes

#ifndef FickianBosanquet_H
#define FickianBosanquet_H


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
                Class FickianBosanquet Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class FickianBosanquet
:
    public MultiSpeciesTransportModel<Thermo>
{
    // Private data

        //- Mean pore diametre
        scalarList dp_;

        // TODO: Should have access to the porosityModel list of the region!
        // Porosity
        scalarList eps_;

        // Tortuosity
        scalarList tau_;

        //- Mass diffusion coefficients in the mixture Fick law - D_iMixture
        PtrList<volScalarField> Dm_;

        //- Effective mass diffusion coefficients in the mixture
        PtrList<volScalarField> DmK_;

        //- Knudsen diffusion coefficients
        PtrList<volScalarField> DK_;  

        //- Molar fractions [-]
        const PtrList<volScalarField>& X_;

        //- Binary diffusivity models 
        PtrList<binaryDiffusivityModel<Thermo>> DijModels_; 

        //- Binary diffusion coefficients
        PtrList<volScalarField> Dij_;  


public:

    //- Runtime type information
    TypeName("FickianBosanquet");

    // Constructors

        //- Construct from a momentum transport model and a thermo model
        FickianBosanquet
        (
            const phaseModel& phase
        );


    //- Destructor
    virtual ~FickianBosanquet();


    // Member Functions

   //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const;

        //- Effective mass diffusion coefficient
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
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FickianBosanquet.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
