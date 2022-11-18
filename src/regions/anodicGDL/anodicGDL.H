/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::anodicGDL.C

SourceFiles
    anodicGDL.C

\*---------------------------------------------------------------------------*/

#ifndef anodicGDL_H
#define anodicGDL_H

#include "regionType.H"
#include "dimensionedScalar.H"
#include "multiRegionSystem.H"
#include "volFields.H"
#include "IOReferencer.H"
#include "HashPtrTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{

/*---------------------------------------------------------------------------*\
                  Class anodicGDL Declaration
\*---------------------------------------------------------------------------*/

class anodicGDL
:
    public regionType
{
        //- name
        word regionName_;

        //- transport properties dictionary 
        IOdictionary transportProperties_;

        //- material properties dictionary 
        IOdictionary materialProperties_;

        //- operating conditions dictionary
        IOdictionary operatingConditions_;

        //- specific heat
        dimensionedScalar cv_;

        //- density
        dimensionedScalar rho_;

        //- thermal conductivity
        autoPtr<volScalarField> k_;

        //- electric conductivity
        autoPtr<volScalarField> sigma_;

        //- diffusion coefficient hydrogen
        dimensionedScalar DH2_;

        //- diffusion coefficient vapor
        dimensionedScalar DV_;

        //- porosity
        dimensionedScalar epsilonP_;

        //- tortuosity
        dimensionedScalar tau_;

        //- opearting pressure
        dimensionedScalar p_;

        //- total gas species concentration
        volScalarField c_;

        //- effective diffusion coefficient hydrogen
        autoPtr<volScalarField> DEffH2_;

        //- effective diffusion coefficient vapor
        autoPtr<volScalarField> DEffV_;

        //- heat sources
        volScalarField sT_;

    // Private member functions	

        //- update gas species transport parameters
        void updateGasSpeciesTransportProperties();

        //- update source terms
        void updateSourceTerms();

protected:

        //- temperature
        autoPtr<volScalarField> T_;

        //- electric Potential
        autoPtr<volScalarField> phiE_;	

        //- hydrogen mole fraction
        autoPtr<volScalarField> xH2_;

        //- vapor mole fraction
        autoPtr<volScalarField> xV_;

public:

    TypeName("anodicGDL");

    // Constructors

        anodicGDL
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~anodicGDL();


    // Member Functions

        // Access
        const volScalarField& T() const
        {
            return T_();
        }

        volScalarField& T()
        {
            return T_();
        }

        const volScalarField& phiE() const
        {
            return phiE_();
        }

        volScalarField& phiE()
        {
            return phiE_();
        }

        const volScalarField& xH2() const
        {
            return xH2_();
        }

        volScalarField& xH2()
        {
            return xH2_();
        }

        const volScalarField& xV() const
        {
            return xV_();
        }

        volScalarField& xV()
        {
            return xV_();
        }

        virtual bool active() const
        {
            return true;
        }

        virtual void correct();

        virtual scalar getMinDeltaT();

        virtual void solveRegion();

        virtual void prePredictor();

        virtual void momentumPredictor();

        virtual void pressureCorrector();

        virtual void meshMotionCorrector();

        virtual void setCoupledEqns();

        virtual void postSolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
