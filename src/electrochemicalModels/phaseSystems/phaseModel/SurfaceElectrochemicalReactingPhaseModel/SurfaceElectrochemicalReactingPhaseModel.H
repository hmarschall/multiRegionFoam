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
    Foam::SurfaceElectrochemicalReactingPhaseModel

Description

SourceFiles
    SurfaceElectrochemicalReactingPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef SurfaceElectrochemicalReactingPhaseModel_H
#define SurfaceElectrochemicalReactingPhaseModel_H

#include "phaseModel.H"
#include "rhoReactionThermo.H"
#include "saturationModel.H"

#include "regionType.H"
// #include "activationOverpotentialModel.H"

#include "patchToPatchInterpolation.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class SurfaceElectrochemicalReactingPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class SurfaceElectrochemicalReactingPhaseModel
:
    public BasePhaseModel
{

private:

    //- Activation overpotential variables

        //- Activation overpotential model
        volScalarField eta_;

        //- Jump (eta + nernst)
        volScalarField jump_;

        //- Catalyst layer adjacent regions
        const dictionary& regions_;

        //- Saturation model
        autoPtr<saturationModel> saturation_;

        // phi names for electric field
        // electron: phi
        // ion: phi
        HashTable<word> phiNames_;

        // i names for electric field
        // electron: i
        // ion: i
        HashTable<word> iNames_;

        //- Dissolved on/off
        Switch dissolved_;

        // Properties
        // relax, alpha, j0, Xref, gamma
        scalar relax_;
        scalar alpha_;
        scalar gamma_;

        // Current exchange density value A/m^2
        dimensionedScalar i0Value_;

        // Current exchange density A/m^2
        volScalarField i0_;

        //- Current density
        //- A/m2
        volScalarField i_;

        // //- Create nernst model
        // void createNernst();

        HashTable<scalar> rxnList_;

    //- Nernst variables

        // Nernst potential
        volScalarField nernst_;

        //- delta H
        volScalarField deltaH_;

        //- delta S
        volScalarField deltaS_;

        //- Reference pressure
        dimensionedScalar pRef_;

        // Species dictionary
        const dictionary& species_;

        //- Coefficient
        scalar A_;
        scalar B_;

        // Heat source due to reactions for visualization
        volScalarField Qdot_;

public:

    // Constructors

        SurfaceElectrochemicalReactingPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const label index
        );


    //- Destructor
    virtual ~SurfaceElectrochemicalReactingPhaseModel();


    // Member Functions

        //- Correct electroChemicalReaction rate
        virtual void correct();

        //- Fuel consumption rate matrix, i.e. source term for fuel equation
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;


        //- Fuel consumption rate matrix, i.e. source term for fuel equation
        void correctGasIdmdt(const volScalarField& Y) const;

        //- Return the heat source/sink (as explicit form)
        virtual tmp<volScalarField> Qdot();

                //- Lookup regionType
        const regionType& region(const word& name) const
        {
            return this->fluid().mesh().time().template
                lookupObject<regionType>(name + "Dict");    // Look how the regionTypes are registered!
        }


    //- Electrochemical member functions
        virtual void correctOverpotential();

        virtual void correctNernst();

        //- J/mol
        const volScalarField& deltaH() const
        {
            return deltaH_;
        }

        //- J/(mol K)
        const volScalarField& deltaS() const
        {
            return deltaS_;
        }

        scalar eta(scalar min, scalar max, scalar i0, scalar i) const;



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SurfaceElectrochemicalReactingPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
