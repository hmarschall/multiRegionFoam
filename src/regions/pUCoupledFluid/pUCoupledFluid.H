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
    Foam::pUCoupledFluid.C

SourceFiles
    pUCoupledFluid.C

\*---------------------------------------------------------------------------*/

#ifndef pUCoupledFluid_H 
#define pUCoupledFluid_H 

#include "regionType.H"
#include "dimensionedScalar.H"
#include "multiRegionSystem.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOReferencer.H"
#include "HashPtrTable.H"
#include "MRFZones.H"
#include "fvCFD.H"
#include "fvBlockMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{

/*---------------------------------------------------------------------------*\
                           Class pUCoupledFluid Declaration
\*---------------------------------------------------------------------------*/

class pUCoupledFluid
:
    public regionType
{       

        //- name
        word regionName_;      

        //- transport
        autoPtr<volVectorField> U_;

        autoPtr<surfaceScalarField> phi_;

        autoPtr<volScalarField> p_;

        //- source term of pressure poisson equation
        autoPtr<surfaceScalarField> presSource_;

        //- Block vector field for velocity (first entry) and pressure (second
        //  entry)
        autoPtr<volVector4Field> Up_;

        //- reverse coefficient field of U for pressure poisson equation
        autoPtr<volScalarField> rAU_;

        //- gradients
        autoPtr<volVectorField> gradp_;
        autoPtr<volTensorField> gradU_;

        //- pcorr
        wordList pcorrTypes_;
        autoPtr<volScalarField> pcorr_;

        //- pressure field bounds
        dimensionedScalar pMin_;
        dimensionedScalar pMax_;

        //- velocity field bounds
        dimensionedScalar UMax_;
        dimensionedScalar smallU_;

        //- transport properties dictionary 
        IOdictionary transportProperties_;

        //- density
        dimensionedScalar rhoFluid_;
        autoPtr<volScalarField> rho_;

        //- viscosity
        dimensionedScalar muFluid_;
        autoPtr<volScalarField> mu_;

        //- closed volume?
        Switch closedVolume_;

        //- has space patch?
        Switch hasSpacePatch_;

        //- Reference pressure cell
        label pRefCell_;

        //- Reference pressure value
        scalar pRefValue_;

        //- List of lable indicating if proc holds pRefCell
        labelList whichProcHasRef_;
        
        //- block solver residual
        // BlockSolverPerformance<vector4> residual_;
        // scalar maxResidual;
        // scalar convergenceCriterion;

        // //- Current corrector loop index
        // label corr_;

        //- Moving reference frames
        MRFZones mrfZones_;
        label myTimeIndex_;

        //- Adjustable time-step switch
        const Switch adjustTimeStep_;

        //- For adjustable time-step, this is the maximum Courant number
        const scalar maxCo_;

        //- For adjustable time-step, this is the maximum time-step
        const scalar maxDeltaT_;

        //- Indicator to correct phi after first solution step
        bool firstUpdate_;

public:

    TypeName("pUCoupledFluid");

    // Constructors

        pUCoupledFluid
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~pUCoupledFluid();


    // Member Functions

        // Access
        const volVectorField& U() const
        {
            return U_();
        }

        const volScalarField& p() const
        {
            return p_();
        }        

        volVectorField& U()
        {
            return U_();
        }

        volScalarField& p()
        {
            return p_();
        }

        // Access by name
        template<class T>
        const T& getObject(const word& name)
        {
            return
            (
                mesh().thisDb().lookupObject<T>(name)
            );
        }

        template<class T>
        T& getObject(const word& name)
        {
            return
            (
                const_cast<const T&>(mesh().thisDb().lookupObject<T>(name))
            );
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
