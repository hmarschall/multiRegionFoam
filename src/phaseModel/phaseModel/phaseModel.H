/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::phaseModel

SourceFiles
    phaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseModel_H
#define phaseModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
#include "basicRhoThermo.H"
#include "runTimeSelectionTables.H"

#include "porosityModelList.H"
#include "regionType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class phaseModel
:
    public volScalarField
{
    // Private data

        //- Mesh of the phase model
        const fvMesh& mesh_;

        //- Dictionary of the phase model
        const dictionary& dict_;

        //- Porous zone list
        // TODO: Should not be here!
        // -> Should be placed in the phaseSystem, because the other phase should not define it
        // too
        autoPtr<porosityModelList> porosityModels_;

        //- Name of the phaseModel
        const word name_;

public:

    //- Runtime type information
    ClassName("phaseModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );


    // Constructors

        phaseModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Return clone
        autoPtr<phaseModel> clone() const;


    // Selectors

        static autoPtr<phaseModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~phaseModel();


    // Member Functions

        //- Return the mesh of this phase
        const fvMesh& mesh() const;

        //- Return the dictionary of this phase
        const dictionary& dict() const;

        //- Lookup regionType
        const regionType& region(const word& name) const
        {
            return mesh().time().
                lookupObject<regionType>(name);
        }

        // Return the name of the phase
        const word& name() const;

        //- Return the name of the phase for use as the keyword in PtrDictionary
        const word& keyword() const;

        //- Correct the phase properties
        virtual void correct();

        //- Correct the kinematics
        virtual void correctKinematics();

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Correct the turbulence
        virtual void correctTurbulence();

        //- Correct the energy transport
        virtual void correctEnergyTransport();

        virtual Foam::porosityModelList& porosityModels()
        {
            return porosityModels_();
        }

        // Compressibility (variable density)

            //- Return true if the phase is compressible otherwise false
            virtual bool compressible() const = 0;

            //- Return the phase dilatation rate (d(alpha)/dt + div(alpha*phi))
            virtual tmp<volScalarField> divU() const = 0;

            //- Set the phase dilatation rate (d(alpha)/dt + div(alpha*phi))
            virtual void divU(tmp<volScalarField> divU) = 0;


        // // Thermo

            //- Return whether the phase is isothermal
            virtual bool isothermal() const = 0;

            //- Return the enthalpy equation
            virtual tmp<fvScalarMatrix> heEqn() = 0;


            //- Return all of the heat source/sink (as explicit form)
            virtual tmp<volScalarField> heQdot() = 0;

            // //- Return the heat source/sink (as explicit form)
            // virtual tmp<volScalarField> Qdot() const = 0;

            //- Return the thermophysical model
            virtual const basicRhoThermo& thermo() const = 0;

            //- Access the thermophysical model
            virtual basicRhoThermo& thermoRef() = 0;

            //- Return the density field
            virtual tmp<volScalarField> rho() const = 0;


        // Species

            //- Return whether the phase is pure (i.e., not multi-component)
            virtual bool pure() const = 0;

            //- Return the species fraction equation
            virtual tmp<fvScalarMatrix> YiEqn(volScalarField& Yi) = 0;

            //- Return the source/sink term for species
            virtual const PtrList<volScalarField>& iDmdt() const = 0;

            //- Return the source/sink term for species by name
            virtual const volScalarField& iDmdt(const word& name) const = 0;

            //- Access the source/sink term for species
            virtual PtrList<volScalarField>& iDmdtRef() = 0;

            //- Return the total mass transfer rate
            virtual tmp<volScalarField> dmdt() const = 0;

            //- Return the species mass fractions
            virtual const PtrList<volScalarField>& Y() const = 0;

            //- Return a species mass fraction by name
            virtual const volScalarField& Y(const word& name) const = 0;

            //- Access the species mass fractions
            virtual PtrList<volScalarField>& YRef() = 0;

            //- Return the active species mass fractions
            virtual const UPtrList<volScalarField>& YActive() const = 0;

            //- Access the active species mass fractions
            virtual UPtrList<volScalarField>& YActiveRef() = 0;

            //- Return the species mass fractions
            virtual const PtrList<volScalarField>& X() const = 0;

            //- Return a species mass fraction by name
            virtual const volScalarField& X(const word& name) const = 0;

            //- Access the species mass fractions
            virtual PtrList<volScalarField>& XRef() = 0;

            //- Return the active species mass fractions
            virtual const UPtrList<volScalarField>& XActive() const = 0;

            //- Access the active species mass fractions
            virtual UPtrList<volScalarField>& XActiveRef() = 0;


        // // Momentum

            //- Return whether the phase is stationary
            virtual bool stationary() const = 0;

            //- Return the momentum equation
            virtual tmp<fvVectorMatrix> UEqn() = 0;

            //- Return the momentum equation for the face-based algorithm
            virtual tmp<fvVectorMatrix> UfEqn() = 0;

            //- Return the velocity
            virtual tmp<volVectorField> U() const = 0;

            //- Access the velocity
            virtual volVectorField& URef() = 0;

            //- Return the volumetric flux
            virtual tmp<surfaceScalarField> phi() const = 0;

            //- Access the volumetric flux
            virtual surfaceScalarField& phiRef() = 0;

            //- Return the volumetric flux of the phase
            virtual tmp<surfaceScalarField> alphaPhi() const = 0;

            //- Access the volumetric flux of the phase
            virtual surfaceScalarField& alphaPhiRef() = 0;

            //- Return the mass flux of the phase
            virtual tmp<surfaceScalarField> alphaRhoPhi() const = 0;

            //- Access the mass flux of the phase
            virtual surfaceScalarField& alphaRhoPhiRef() = 0;

            //- Return the continuity error
            virtual tmp<volScalarField> continuityError() const = 0;

            //- Return the continuity error due to the flow field
            virtual tmp<volScalarField> continuityErrorFlow() const = 0;

            //- Return the continuity error due to any sources
            virtual tmp<volScalarField> continuityErrorSources() const = 0;

            //- Return the phase kinetic energy
            virtual tmp<volScalarField> K() const = 0;


        // // Transport

            //- Return the laminar dynamic viscosity
            virtual tmp<volScalarField> mu() const = 0;

            // //- Return the laminar dynamic viscosity on a patch
            // virtual tmp<scalarField> mu(const label patchi) const = 0;

            //- Return the laminar kinematic viscosity
            virtual tmp<volScalarField> nu() const = 0;

            // //- Return the laminar kinematic viscosity on a patch
            // virtual tmp<scalarField> nu(const label patchi) const = 0;

            //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
            virtual tmp<volScalarField> alpha() const = 0;

            // //- Thermal diffusivity for enthalpy of mixture for patch [kg/m/s]
            // virtual tmp<scalarField> alpha(const label patchi) const = 0;

            //- Thermal diffusivity for temperature of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappa() const = 0;

            // //- Thermal diffusivity for temperature of mixture
            // //  for patch [J/m/s/K]
            // virtual tmp<scalarField> kappa(const label patchi) const = 0;

        //     //- Thermal diffusivity for energy of mixture [kg/m/s]
        //     virtual tmp<volScalarField> alphahe() const = 0;

        //     //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
        //     virtual tmp<scalarField> alphahe(const label patchi) const = 0;


        // Electrochemistry

            // //- Fuel consumption rate matrix, i.e. source term for fuel equation
            // virtual void correctGasIdmdt(const volScalarField& Y) const = 0;

            virtual void correctOverpotential() = 0;

            virtual void correctNernst() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
