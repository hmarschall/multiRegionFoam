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
    Foam::catalystLayer.C

SourceFiles
    catalystLayer.C

\*---------------------------------------------------------------------------*/

#ifndef catalystLayer_H
#define catalystLayer_H

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
                  Class catalystLayer Declaration
\*---------------------------------------------------------------------------*/

class catalystLayer
:
    public regionType
{
        //- mesh
        const fvMesh& mesh_;

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

	//- permeability
	dimensionedScalar K0_;

        //- thermal conductivity
        autoPtr<volScalarField> k_;

	//- electric conductivity
	autoPtr<volScalarField> sigma_;
	
	//- protonic conductivity
	autoPtr<volScalarField> kappa_;

	//- diffusion coefficient oxygen
	dimensionedScalar DO2_;

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

	//- irreducible saturation
	dimensionedScalar sIm_;

	//- liquid gas exchange area
	dimensionedScalar aLG_;

	//- pt Loading
	dimensionedScalar a_;

	//- symmetry factor
	dimensionedScalar beta_;

	//- reaction entropy
	dimensionedScalar deltaS_;

	//- reaction enthalpy
	dimensionedScalar deltaH_;

	//- reaction activation energy
	dimensionedScalar ER_;

	//- pre factor exchange current density
	dimensionedScalar j0_;

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

	//- effective diffusion coefficient oxygen
	autoPtr<volScalarField> DEffO2_;

	//- effective diffusion coefficient vapor
	autoPtr<volScalarField> DEffV_;

	//- reduced liquid water saturation
	volScalarField sRed_;

	//- saturation vapor fraction
	volScalarField xVSat_;

	//- permeability
	autoPtr<volScalarField> K_;

	//- dynamic viscosity of liquid water
	volScalarField mu_;

	//- capillary pressure derived after saturation
	volScalarField dpCds_;

	//- evaporation/condensation rate
	volScalarField gamma_;

	//- overpotential
	volScalarField eta_;

	//- volumetric exchange current density
	volScalarField j_;

	//- heat sources
	volScalarField sT_;

	//- mass source water in ionomer
	volScalarField sLambda_;

	//- mass source vapor
	volScalarField sV_;

	//- mass source liquid water
	volScalarField ss_;

	

protected:

        //- temperature
//        volScalarField T_;

        autoPtr<volScalarField> T_;

	//- electric Potential
//	  volScalarField phiE_;

	autoPtr<volScalarField> phiE_;	

	//- electrolyte Potential
//	  volScalarField phiP_;

	autoPtr<volScalarField> phiP_;

	//- water content in ionomer
//	  volScalarField lambda_;

	autoPtr<volScalarField> lambda_;

	//- oxygen mole fraction
//	  volScalarField xO2_;

	autoPtr<volScalarField> xO2_;

	//- vapor mole fraction
//	  volScalarField xv_;

	autoPtr<volScalarField> xV_;

	//- saturation
//	  volScalarField s_;

	autoPtr<volScalarField> s_;

public:

    TypeName("catalystLayer");

    // Constructors

        catalystLayer
        (
            const fvMesh& mesh,
            const word& regionName
        );

    //- Destructor
    virtual ~catalystLayer();


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

	const volScalarField& xO2() const
        {
            return xO2_();
        }

        volScalarField& xO2()
        {
            return xO2_();
        }

	const volScalarField& xV() const
	{
	    return xV_();
	}

	volScalarField& xV()
	{
	    return xV_();
	}

	const volScalarField& s() const
	{
	    return s_();
	}

	volScalarField& s()
	{
	    return s_();
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
