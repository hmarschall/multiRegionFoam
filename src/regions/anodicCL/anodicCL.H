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
    Foam::anodicCL.C

SourceFiles
    anodicCL.C

\*---------------------------------------------------------------------------*/

#ifndef anodicCL_H
#define anodicCL_H

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
                  Class anodicCL Declaration
\*---------------------------------------------------------------------------*/

class anodicCL
:
    public regionType
{
        //- name
        word regionName_;

        //- transport properties dictionary 
        IOdictionary transportProperties_;

        //- material properties dictionary 
        IOdictionary materialProperties_;

        //- electrochemical properties dictionary
        IOdictionary electrochemicalProperties_;

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

        //- protonic conductivity
        autoPtr<volScalarField> kappa_;

        //- diffusion coefficient hydrogen
        dimensionedScalar DH2_;

        // diffusion coefficient vapor
        dimensionedScalar DV_;

        //- porosity
        dimensionedScalar epsilonP_;

        //- tortuosity
        dimensionedScalar tau_;

        //- ionomer content
        dimensionedScalar epsilonI_;

        //- layer thickness
        dimensionedScalar d_;

        //- molecular weight water;
        dimensionedScalar MW_;

        //- molar volume of liquid water
        dimensionedScalar VW_;

        //- molar volume of membrane
        dimensionedScalar VM_;

        //- pt Loading
        dimensionedScalar a_;

        //- symmetry factor
        dimensionedScalar beta_;

        //- reaction entropy
        dimensionedScalar deltaS_;

        //- reaction activation energy
        dimensionedScalar ER_;

        //- pre factor exchange current density
        dimensionedScalar jStern_;

        //- opearting pressure
        dimensionedScalar p_;

        //- relative humidity
        dimensionedScalar RH_;

        //- diffusion coefficient water in ionomer
        autoPtr<volScalarField> DLambda_;

        //- volume fraction of water in membrane
        volScalarField f_;

        //- volume fraction of water in membrane for protonic conductivity
        volScalarField fCond_;

        //- electro-osmotic drag coefficient
        volScalarField xi_;

        //- equilibrium water content in ionomer
        volScalarField lambdaEq_;

        //- ab-/desorption rate
        volScalarField kSorp_;

        //- total gas species concentration
        volScalarField c_;

        //- effective diffusion coefficient hydrogen
        autoPtr<volScalarField> DEffH2_;

        //- effective diffusion coefficient vapor
        autoPtr<volScalarField> DEffV_;

        //- saturation vapor fraction
        volScalarField xVSat_;

        //- reaction exchange current density
        volScalarField j0_;

        //- Galvani potential difference
        volScalarField deltaPhi_;

        //- reversible potential difference
        volScalarField deltaPhi0_;

        //- overpotential
        volScalarField eta_;

        //- volumetric current density source
        volScalarField j_;

        //- heat sources
        volScalarField sT_;

        //- mass source water in ionomer
        volScalarField sLambda_;

        //- mass source vapor
        volScalarField sV_;

    // Private member functions

        //- update ionomer properties
	void updateIonomerProperties();

        //- update gas species transport parameters
        void updateGasSpeciesTransportProperties();

        //- update ab-/desorption
        void updateAbsorptionDesorption();

        //- update electrochemistry
        void updateElectrochemistry();

        //- update source terms
        void updateSourceTerms();

protected:

        //- temperature
        autoPtr<volScalarField> T_;

        //- electric Potential
        autoPtr<volScalarField> phiE_;	

        //- electrolyte Potential
        autoPtr<volScalarField> phiP_;

        //- water content in ionomer
        autoPtr<volScalarField> lambda_;

        //- hydrogen mole fraction
        autoPtr<volScalarField> xH2_;

	//- vapor mole fraction
        autoPtr<volScalarField> xV_;

public:

    TypeName("anodicCL");

    // Constructors

        anodicCL
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~anodicCL();


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

        const volScalarField& phiP() const
        {
            return phiP_();
        }

        volScalarField& phiP()
        {
            return phiP_();
        }

        const volScalarField& lambda() const
        {
            return lambda_();
        }

        volScalarField& lambda()
        {
            return lambda_();
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

        virtual void setRDeltaT();

        virtual void solveRegion();

        virtual void setCoupledEqns();

        virtual void updateFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
