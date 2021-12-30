/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    interfaceCoupledPressureFlux

Author
    Holger Marschall (marschall@csi.tu-darmstadt.de)

Description
    Neumann boundary condition for pressure within the ALE interface tracking.
    Normal pressure flux jump is calculated from the neighbour patch flux.

    Example usage:
        interfacePatch
            {
                type                    interfaceCoupledPressureFlux;
                k                       "rho";
                neighbourRegionName     liquid;
                neighbourPatchName      interfacePatchShadow;
                value                   uniform 0;
            }


SourceFiles
    interfaceCoupledPressureFlux.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceCoupledPressureFlux_H
#define interfaceCoupledPressureFlux_H

#include "fixedGradientFvPatchFields.H"
#include "autoPtr.H"

#include "patchCoupleManager.H"
#include "regionInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class interfaceCoupledPressureFlux Declaration
\*---------------------------------------------------------------------------*/

class interfaceCoupledPressureFlux
:
    public fixedGradientFvPatchScalarField, 
    public patchCoupleManager    
{
    // Private data

        //- Diffusion coefficient
		word kName_;
		
	    //- Neighbour names
        word neighbourRegionName_;

        word neighbourPatchName_;

        //- Fluid properties
        dimensionedScalar muFluidA_;

        dimensionedScalar muFluidB_;

        dimensionedScalar rhoFluidA_;

        dimensionedScalar rhoFluidB_;

        //- Gravitational accelaration
        dimensionedVector g_;

        //- Calculate pressure gradient jump
        tmp<scalarField> nGradJump() const;
        

public:

    //- Runtime type information
    TypeName("interfaceCoupledPressureFlux");

    // Constructors

        //- Construct from patch and internal field
        interfaceCoupledPressureFlux
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        interfaceCoupledPressureFlux
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given interfaceCoupledPressureFlux
        //  onto a new patch
        interfaceCoupledPressureFlux
        (
            const interfaceCoupledPressureFlux&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
            	new interfaceCoupledPressureFlux(*this)
            );
        }

        //- Construct as copy setting internal field reference
        interfaceCoupledPressureFlux
        (
            const interfaceCoupledPressureFlux&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new interfaceCoupledPressureFlux(*this, iF)
            );
        }


    // Member functions

        //- Update the patch field coefficients
        virtual void updateCoeffs();

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::blocking
        );

        //- Return the maximum normalized coupled patch residual
	    scalarField residual() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
