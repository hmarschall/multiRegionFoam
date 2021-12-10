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
    Foam::fluidFluid.C

SourceFiles
    fluidFluid.C

\*---------------------------------------------------------------------------*/

#ifndef fluidFluid_H
#define fluidFluid_H

#include "regionInterface.H"
#include "dimensionedScalar.H"
#include "multiRegionSystem.H"
#include "volFields.H"
#include "IOReferencer.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionInterfaces
{

/*---------------------------------------------------------------------------*\
                           Class fluidFluid Declaration
\*---------------------------------------------------------------------------*/

class fluidFluid
:
    public regionInterface
{
        IOdictionary transportPropertiesA_;

        IOdictionary transportPropertiesB_;

        IOdictionary gravitationalProperties_;

        //const Time& runTime_;
        //const fvPatch& patchA_;
        //const fvPatch& patchB_;

        //- bulk velocity field        
        const volVectorField& U_;

        //- bulk flux field
        const surfaceScalarField& phi_;
        
        //- phase densities
        dimensionedScalar rhoA_;

        dimensionedScalar rhoB_;

        //- phase viscosities
        dimensionedScalar muA_;

        dimensionedScalar muB_;
        
        //- surface tension coefficient
        dimensionedScalar sigma0_;

        //- gravity
        const dimensionedVector g_; 
              
public:

    TypeName("fluidFluid");

    // Constructors

        fluidFluid
        (
            const Time& runTime,
            const fvPatch& patchA,
            const fvPatch& patchB
        );

    // Destructor
    virtual ~fluidFluid()
    {}
    
      
    // Access  
    
        //- Return reference to velocity
        const volVectorField& U() const
        {
            return U_;
        }
 
    // Member Functions	
    
		//- Return constant reference to velocity field on patch
        const vectorField& Up();
        
        //- Update surface velocity field
        void updateUs();

        //- Update surface velocity b.c.
        void correctUsBoundaryConditions();

        //- Update surface flux
        void updateSurfaceFlux();

        //- Update interfacial curvature field
        void updateK();

        //- Calculate interface curvature at axis
        void calcCurvatureAxis
        (
            areaScalarField& K
        );
    
    // virtual functions 
        
        //- Create surface velocity field
        virtual void makeUs() const;
        
        //- Create interfacial curvature                 
        virtual void makeK() const;

        //- Create surface flux
        virtual void makePhis() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionInterfaces
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
