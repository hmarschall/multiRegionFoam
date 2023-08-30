/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::basicSpecieMixture

Group
    grpReactionThermophysicalMixtures

Description
    Specialization of basicMultiComponentMixture for a mixture consisting
    of a number for molecular species.

SourceFiles
    basicSpecieMixture.C

\*---------------------------------------------------------------------------*/

#ifndef basicSpecieMixture_H
#define basicSpecieMixture_H

// #include "basicMultiComponentMixture.H"
#include "volFields.H"
#include "PtrList.H"
#include "speciesTable.H"
#include "basicMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class basicSpecieMixture Declaration
\*---------------------------------------------------------------------------*/

class basicSpecieMixture
:
    public basicMixture
{


protected:

    // Protected data

        // TODO: I am not allowed to define this here as well!! -> is in Ofdev but here it seems to be somewhere else
        // //- Phase-name
        // const word& phaseName_;

        //- Table of specie names
        speciesTable species_;

        //- The name of the default specie
        //  The mass fraction of which is derived from the other species
        //  rather than solved.
        //  Also used as the carrier specie in multicomponent diffusion
        word defaultSpecie_;

        //- The index of the default specie
        label defaultSpecieIndex_;

        //- List of specie active flags
        mutable List<bool> active_;

        //- Species mass fractions [-]
        PtrList<volScalarField> Y_;

        // //- Species mole fractions [-]
        // PtrList<volScalarField> X_;

        // //- Species molar concentrations [mol/m^3]
        // PtrList<volScalarField> C_;

public:

    //- Run time type information
    TypeName("basicSpecieMixture");

    //- The base class of the mixture
    typedef basicSpecieMixture basicMixtureType;


    // Constructors

        //- Construct from dictionary, species names, mesh and phase name
        basicSpecieMixture
        (
            const dictionary& thermoDict,
            const wordList& specieNames,
            const fvMesh& mesh,
            const word& phaseName
        );


    //- Destructor
    virtual ~basicSpecieMixture() = default;


    // Member functions

    //- Return the table of species
        const speciesTable& species() const
        {
            return species_;
        }

        //- Does the mixture include this specie?
        inline bool contains(const word& specieName) const;

        //- Return the index of the default specie
        inline label defaultSpecie() const;

        //- Return true for active species
        inline bool active(label speciei) const;

        //- Return the bool list of active species
        inline const List<bool>& active() const;

        //- Set speciei active
        inline void setActive(label speciei);

        //- Set speciei inactive
        inline void setInactive(label speciei);

        //- Return the mass-fraction fields
        inline PtrList<volScalarField>& Y();

        //- Return the const mass-fraction fields
        inline const PtrList<volScalarField>& Y() const;

        //- Return the mass-fraction field for a specie given by index
        inline volScalarField& Y(const label i);

        //- Return the const mass-fraction field for a specie given by index
        inline const volScalarField& Y(const label i) const;

        //- Return the mass-fraction field for a specie given by name
        inline volScalarField& Y(const word& specieName);

        //- Return the const mass-fraction field for a specie given by name
        inline const volScalarField& Y(const word& specieName) const;

        // //- Return the mole fraction fields
        // inline PtrList<volScalarField>& X();

        // //- Return the const mole fraction fields
        // inline const PtrList<volScalarField>& X() const;

        // //- Return the mole fraction field for a specie given by index
        // inline volScalarField& X(const label i);

        // //- Return the const mole fraction field for a specie given by index
        // inline const volScalarField& X(const label i) const;

        // //- Return the mole fraction field for a specie given by name
        // inline volScalarField& X(const word& specieName);

        // //- Return the const mole fraction field for a specie given by name
        // inline const volScalarField& X(const word& specieName) const;

        // //- Return the molar concentration fields
        // inline PtrList<volScalarField>& C();

        // //- Return the const molar concentration fields
        // inline const PtrList<volScalarField>& C() const;

        // //- Return the molar concentration field for a specie given by index
        // inline volScalarField& C(const label i);

        // //- Return the const molar concentration field for a specie given by index
        // inline const volScalarField& C(const label i) const;

        // //- Return the molar concentration field for a specie given by name
        // inline volScalarField& C(const word& specieName);

        // //- Return the const molar concentration field for a specie given by name
        // inline const volScalarField& C(const word& specieName) const;

        // //- Update molar fields (mole fractions and molar concentrations)
        // inline void updateMolarFields();

        //- Return the specie index of the given mass-fraction field
        inline label index(const volScalarField& Yi) const;

          // Specie properties

            //- Molecular weight of the given specie [kg/kmol]
            virtual scalar Wi(const label speciei) const = 0;

            //- Molar volume of the given specie [m^3/mol]
            virtual scalar Vi(const label speciei) const = 0;

            //- Valence number of the given specie [kg/kmol]
            virtual scalar Zi(const label speciei) const = 0;

            // //- Enthalpy of formation [J/kg]
            // virtual scalar Hf(const label speciei) const = 0;


        // Specie thermo properties

            //- Density [kg/m^3]
            virtual scalar rho
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual scalar Cp
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual scalar HE
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<scalarField> HE
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<volScalarField> HE
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual scalar Hs
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual tmp<scalarField> Hs
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual tmp<volScalarField> Hs
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Entropy [J/(kg K)]
            virtual scalar S
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Entropy [J/(kg K)]
            virtual tmp<scalarField> S
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T
            ) const = 0;

            //- Entropy [J/(kg K)]
            virtual tmp<volScalarField> S
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual scalar Ha
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<scalarField> Ha
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual tmp<volScalarField> Ha
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;


        // Specie transport properties

            //- Dynamic viscosity [kg/m/s]
            virtual scalar mu
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Dynamic viscosity [kg/m/s]
            virtual tmp<volScalarField> mu
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Thermal conductivity [W/m/K]
            virtual scalar kappa
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Thermal conductivity [W/m/K]
            virtual tmp<volScalarField> kappa
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicSpecieMixtureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
