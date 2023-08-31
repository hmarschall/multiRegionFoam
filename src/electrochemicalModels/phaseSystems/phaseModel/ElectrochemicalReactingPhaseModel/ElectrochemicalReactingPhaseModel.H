/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Fork:     foam-extend
    \\  /    A nd           | Version:  4.1                                 
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
    Foam::SurfaceElectrochemicalReactingPhaseModel.C

SourceFiles
    SurfaceElectrochemicalReactingPhaseModel.C

Author
    Steffen Hess(s.hess@fz-juelich.de, Affiliation A) 
    All rights reserved.

Contact
    Steffen Hess(s.hess@fz-juelich.de) 

Affiliations
    Affiliation A)
    Forschungszentrum Jülich IEK-14

Acknowledgement
    Funded by
    AIDAS – AI, Data Analytics and Scalable Simulation

\*---------------------------------------------------------------------------*/

#ifndef ElectrochemicalReactingPhaseModel_H
#define ElectrochemicalReactingPhaseModel_H

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
                  Class ElectrochemicalReactingPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class ElectrochemicalReactingPhaseModel
:
    public BasePhaseModel
{

private:

        //- Activation overpotential model
        volScalarField eta_;

        //- Acitve zone name
        word zoneName_;

        //- Catalyst layer adjacent regions
        const dictionary& regions_;

        //- Saturation model
        autoPtr<saturationModel> saturation_;

        // phi names for electric field
        // electron: phi
        // ion: phi
        HashTable<word> phiNames_;

        //- Dissolved on/off
        Switch dissolved_;

        // Properties
        // relax, alpha, j0, Xref, gamma
        scalar relax_;
        scalar alpha_;
        scalar gamma_;

        dimensionedScalar j0Value_;

        //- Current density
        //- A/m3
        volScalarField j_;

        //- Current exchange density
        //- A/m3
        volScalarField j0_;

        // //- Create nernst model
        // void createNernst();

        //- Activation energy
        dimensionedScalar EA_;

        //- Reference Temperature
        dimensionedScalar Tref_;

        HashTable<scalar> rxnList_;

    //- Nernst variables

        // Nernst potential
        volScalarField nernst_;

        //- delta H
        volScalarField deltaH_;

        //- delta S
        volScalarField deltaS_;

        //- residual Y
        scalar residualY_;

        //- Reference pressure
        dimensionedScalar pRef_;

        // Species dictionary
        // e.g.
        // species
        // {
        //     O2
        //     {
        //         ksi     1.0;            // concentration exponent
        //         cRef    15.5;           // reference concentration
        //     }
        // }
        const dictionary& species_;



public:

    // Constructors

        ElectrochemicalReactingPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const label index
        );


    //- Destructor
    virtual ~ElectrochemicalReactingPhaseModel();


    // Member Functions

        //- Correct electroChemicalReaction rate
        virtual void correct();

        //- Fuel consumption rate matrix, i.e. source term for fuel equation
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Return the heat source/sink (as explicit form)
        virtual tmp<volScalarField> Qdot() const;

                //- Lookup regionType
        const regionType& region(const word& name) const
        {
            return this->fluid().mesh().time().template
                lookupObject<regionType>(name + "Dict");    // Look how the regionTypes are registered!
        }


    //- Electrochemical member functions
        virtual void correctOverpotential();

        virtual void correctNernst();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ElectrochemicalReactingPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
