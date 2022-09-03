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
    Foam::icoFluid.C

SourceFiles
    icoFluid.C

\*---------------------------------------------------------------------------*/

#ifndef icoFluid_H 
#define icoFluid_H 

#include "regionType.H"
#include "MRFZones.H"
#include "dimensionedScalar.H"
#include "multiRegionSystem.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOReferencer.H"
#include "HashPtrTable.H"
#include "pimpleControl.H"
#include "solutionControl.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{

/*---------------------------------------------------------------------------*\
                           Class icoFluid Declaration
\*---------------------------------------------------------------------------*/

class icoFluid
:
    public regionType
{       

        //- name
        word regionName_;

        //- transport properties dictionary 
        IOdictionary transportProperties_;

        //- PIMPLE control
        pimpleControl pimple_;

        //- density
        dimensionedScalar rhoFluid_;

        autoPtr<volScalarField> rho_;

        //- viscosity
        dimensionedScalar muFluid_;

        autoPtr<volScalarField> mu_;           

        //- transport
        word velocityName_;
        autoPtr<volVectorField> U_;

        autoPtr<surfaceScalarField> phi_;

        autoPtr<surfaceScalarField> phiHbyA_;

        autoPtr<volScalarField> p_;

        //- fields for poison equation for the pressure
        autoPtr<volScalarField> AU_;

        autoPtr<volVectorField> HU_;

        //- gradients
        autoPtr<volVectorField> gradp_;

        autoPtr<volTensorField> gradU_;

        //- pcorr
        wordList pcorrTypes_;

        autoPtr<volScalarField> pcorr_;

        //- Under-relaxation factor
        const scalar UUrf_;

        //- closed volume?
        Switch closedVolume_;

        //- has space patch?
        Switch hasSpacePatch_;

        //- Reference pressure cell
        label pRefCell_;

        //- Reference pressure value
        scalar pRefValue_;

        //- residual for inner correctors
        scalar innerResidual_;

        //- residual for pressure outer correctors
        scalar residualPressure_;

        //- Current corrector loop index
        label corr_;

        //- Current non-orthogonal corrector loop index
        label corrNonOrtho_;

        //- Moving reference frames
        MRFZones mrfZones_;
        label myTimeIndex_;

        //- Adjustable time-step switch
        const Switch adjustTimeStep_;

        //- For adjustable time-step, this is the maximum Courant number
        const scalar maxCo_;

        //- For adjustable time-step, this is the maximum time-step
        const scalar maxDeltaT_;

        //- Cached momentum matrix
        //  shared between the momentum predictor and pressure corrector
        tmp<fvVectorMatrix> tUEqn;

public:

    TypeName("icoFluid");

    // Constructors

        icoFluid
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~icoFluid();


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

        virtual void setCoupledEqns();

        virtual void solveRegion();

        virtual void prePredictor();

        virtual void momentumPredictor();

        virtual void pressureCorrector();
        
        virtual void postSolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //