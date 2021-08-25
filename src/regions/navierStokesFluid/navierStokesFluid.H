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
        const fvMesh& mesh_;
        
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
        
        //- surface Properties dictionary
        IOdictionary surfaceProperties_;

        //- thermal conductivity
        dimensionedScalar k_;

        //- specific heat
        dimensionedScalar cp_;

        //- density
        dimensionedScalar rhoFluid_;
        
        volScalarField rho_;
        
        //- temperature
        volScalarField T_;
        
        //- viscosity
        dimensionedScalar muFluid_;
        
        volScalarField mu_;
        
        //- surface tension
        dimensionedScalar cleanSurfaceTension_;
        
        areaScalarField surfaceTension_;              
        
        //- fields for poison equation for the pressure
        volScalarField AU_;
        
        volVectorField HU_;
        
        //- gradients
        volVectorField gradp_;
        
        volTensorField gradU_;
        
        //- pcorr
        wordList pcorrTypes_;
        
        volScalarField pcorr_;
        
public:

    TypeName("navierStokesFluid");

    // Constructors

        navierStokesFluid
        (
            const fvMesh& mesh,          
            const word& regionName
        );

    //- Destructor
    virtual ~navierStokesFluid();


    // Member Functions

        // Access
        const volScalarField& T() const
        {
            return T_;
        }

        volScalarField& T()
        {
            return T_;
        }

        virtual bool active() const
        {
            return true;
        }

        virtual void correct();

        virtual void setRDeltaT();

        virtual void solveRegion();

        virtual void solveCoupledPartitioned();

        // virtual tmp<fvScalarMatrix> coupledFvScalarMatrix() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
