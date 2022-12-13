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
    Foam::ionomer.C

SourceFiles
    ionomer.C

\*---------------------------------------------------------------------------*/

#ifndef ionomer_H
#define ionomer_H

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
                  	Class ionomer Declaration
\*---------------------------------------------------------------------------*/

class ionomer
:
    public regionType
{
        //- name
        word regionName_;

        //- transport properties dictionary 
        IOdictionary transportProperties_;

        //- material properties dictionary 
        IOdictionary materialProperties_;

        //- operating conditions
        IOdictionary operatingConditions_;

        //- specific heat
        dimensionedScalar cv_;

        //- density
        dimensionedScalar rho_;	

        //- thermal conductivity
        autoPtr<volScalarField> k_;

        //- protonic conductivity
        autoPtr<volScalarField> kappa_;

        //- molar volume of liquid water
        dimensionedScalar VW_;

        //- molar volume of membrane
        dimensionedScalar VM_;

        //- relative humidity 
        dimensionedScalar RH_;

        //- diffusion coefficient of water in ionomer
        autoPtr<volScalarField> DLambda_;

        //-  volume fraction of water in membrane
        volScalarField f_;

        //- volume fraction of water in membrane for protonic conductivity
        volScalarField fCond_;

        //- electro-osmotic drag coefficient
        volScalarField xi_;

        //- heat sources
        volScalarField sT_;

        //- water content source - electro-osmotic drag
        volScalarField sLambda_;

        //- cached temperature matrix
        tmp<fvScalarMatrix> TEqn;

        //- cached electrolyte potential matrix
        tmp<fvScalarMatrix> phiPEqn;

        //- cached water content in ionomer matrix
        tmp<fvScalarMatrix> lambdaEqn;

    // Private member functions

        //- update Ionomer Properties
        void updateIonomerProperties();

        //- update Source Terms
        void updateSourceTerms();

protected:

        //- temperature
        autoPtr<volScalarField> T_;

        //- electrolyte potential
        autoPtr<volScalarField> phiP_;	

        //- water content in ionomer
        autoPtr<volScalarField> lambda_;
	

public:

    TypeName("ionomer");

    // Constructors

        ionomer
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~ionomer();


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
