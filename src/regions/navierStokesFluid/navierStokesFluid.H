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
    Foam::navierStokesFluid.C

SourceFiles
    navierStokesFluid.C

\*---------------------------------------------------------------------------*/

#ifndef navierStokesFluid_H 
#define navierStokesFluid_H 

#include "regionType.H"
#include "dimensionedScalar.H"
#include "multiRegionSystem.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOReferencer.H"
#include "HashPtrTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{

/*---------------------------------------------------------------------------*\
                           Class navierStokesFluid Declaration
\*---------------------------------------------------------------------------*/

class navierStokesFluid
:
    public regionType
{
        //- mesh
        const dynamicFvMesh& mesh_;
        
        //- finite area mesh
        faMesh aMesh_;

        //- name
        word regionName_;

        //- transport
        volVectorField U_;

        surfaceScalarField phi_;            

        volScalarField p_;
        
        //- transport properties dictionary 
        IOdictionary transportProperties_;

        //- density
        dimensionedScalar rhoFluid_;
        
        volScalarField rho_;
        
        //- viscosity
        dimensionedScalar muFluid_;
        
        volScalarField mu_;           
        
        //- fields for poison equation for the pressure
        volScalarField AU_;
        
        volVectorField HU_;

        surfaceScalarField phiHbyA_;
        
        //- gradients
        volVectorField gradp_;
        
        volTensorField gradU_;
        
        //- pcorr
        wordList pcorrTypes_;
        
        volScalarField pcorr_;
        
        //- Under-relaxation factor
        const scalar UUrf_;
        
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
      
public:

    TypeName("navierStokesFluid");

    // Constructors

        navierStokesFluid
        (
            const dynamicFvMesh& mesh,          
            const word& regionName
        );

    //- Destructor
    virtual ~navierStokesFluid();


    // Member Functions

        // Access
        const volVectorField& U() const
        {
            return U_;
        }

        const volScalarField& p() const
        {
            return p_;
        }        

        volVectorField& U()
        {
            return U_;
        }

        volScalarField& p()
        {
            return p_;
        }

        // Access by name
        template<class T>
        const T& getObject(const word& name)
        {
            return
            (
                this->thisDb().lookupObject<T>(name)
            );
        }

        template<class T>
        T& getObject(const word& name)
        {
            return
            (
                const_cast<const T&>(this->thisDb().lookupObject<T>(name))
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
