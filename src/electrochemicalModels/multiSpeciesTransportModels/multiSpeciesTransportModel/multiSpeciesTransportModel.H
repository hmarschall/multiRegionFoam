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
    Foam::multiSpeciesTransportModel

Description
    Abstract base class for multispecies mass transport model models.

SourceFiles
    multiSpeciesTransportModelI.H
    multiSpeciesTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef multiSpeciesTransportModel_H
#define multiSpeciesTransportModel_H

#include "porosityModelList.H"
#include "rhoReactionThermo.H"
#include "phaseSystem.H"
#include "phaseModel.H"
#include "runTimeSelectionTables.H"
// #include "diffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                    Class multiSpeciesTransportModel Declaration
\*---------------------------------------------------------------------------*/

class multiSpeciesTransportModel
:
    public IOdictionary
{
    
protected:  
  
    // Protected data

        //- phaseModel
        const phaseModel& phase_;

        //- phaseSystem
        const phaseSystem& phaseSys_;

        //- Porous zones
        const porosityModelList& pZones_ ;

        //- Mesh
        const fvMesh& mesh_;
        


private:
  
    // Private member functions
      
        //- Disallow copy construct
        multiSpeciesTransportModel(const multiSpeciesTransportModel&);

        //- Disallow default bitwise assignment
        void operator=(const multiSpeciesTransportModel&);

public:
 
    //- Runtime type information
    TypeName("multiSpeciesTransportModel");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        multiSpeciesTransportModel,
        dictionary,
        (
            const phaseModel& phase
        ),
        (phase)
    );

    // Constructors	  

        //- Construct from components
        multiSpeciesTransportModel
        (
            const phaseModel& phase
        );

    // Selectors	  
  
        //- Return a reference to the selected multispecies transport model

        static autoPtr<multiSpeciesTransportModel> New
        (
            const phaseModel& phase
        );

    //- Destructor

        virtual ~multiSpeciesTransportModel()
        {}


    // Member Functions

        //- Access function to phase fraction
        const phaseModel& phase() const
        {
            return phase_;
        }

     //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        virtual tmp<volScalarField> alphaEff() const = 0;

        //- Effective thermal turbulent diffusivity of mixture
        //  for patch [kg/m/s]
        virtual tmp<scalarField> alphaEff(const label patchi) const = 0;

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture [W/m/K]
        virtual tmp<volScalarField> kappaEff() const = 0;

        // //- Effective thermal turbulent diffusivity for temperature
        // //  of mixture for patch [W/m/K]
        // virtual tmp<scalarField> kappaEff(const label patchi) const = 0;

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const = 0;

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction for patch [kg/m/s]
        virtual tmp<scalarField> DEff
        (
            const volScalarField& Yi,
            const label patchi
        ) const = 0;

        //- Return the heat flux [W/m^2]
        virtual tmp<surfaceScalarField> q() const = 0;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const = 0;

        //- Return the specie flux for the given specie mass-fraction [kg/m^2/s]
        virtual tmp<surfaceScalarField> j(const volScalarField& Yi) const = 0;

        //- Return the source term for the given specie mass-fraction equation
        virtual tmp<fvScalarMatrix> divj(volScalarField& Yi) const = 0;

        //- Solve the thermophysical transport model equations
        //  and correct the transport coefficients
        virtual void correct() = 0;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
