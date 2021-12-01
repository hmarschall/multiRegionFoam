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
        const Time& runTime_;
        const fvPatch& patchA_;
        const fvPatch& patchB_;
        
        //- fvMesh A
        const fvMesh& meshA_; 
        
        //- fvMesh B
        const fvMesh& meshB_;

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
        dimensionedVector g_; 
        
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
    virtual ~fluidFluid();
        
   
    // Member Functions
      
    // Access  
        //- Return reference to velocity
        const volVectorField& U() const
        {
            return U_;
        }
 
		//- Return reference to fvMesh A
		const fvMesh& meshA() const
		{
			return meshA_;
		}   
		
		//- Return reference to fvMesh B     
		const fvMesh& meshB() const
		{
			return meshB_;
		}     
    
// virtual functions from regionInterface.H
        
        // Update interpolator
        virtual void updateInterpolatorAndGlobalPatches(); 

        // Attach the mesh interface using mesh modifier
        virtual void attach();

        // Detach the mesh interface using mesh modifier
        virtual void detach();  
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionInterfaces
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
