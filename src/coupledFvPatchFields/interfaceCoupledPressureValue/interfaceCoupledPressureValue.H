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
    interfaceCoupledPressureValue

Author
    Holger Marschall (marschall@csi.tu-darmstadt.de)

Description
    Dirichlet boundary condition for pressure within the ALE interface tracking.
    Pressure jump across the interface is accounted for by adding the jump to 
    the neighbour patch value. Normal pressure flux is calculated and divided 
    by the density.

    Example usage:
        interfacePatchShadow
            {
                type                    interfaceCoupledPressureValue;
	            k			            "rho";
	            neighbourRegionName     gas;
                neighbourPatchName      interfacePatch;
                nonOrthCorr             yes;
                secondOrder             yes;
                value                   uniform 0;
            }


SourceFiles
    interfaceCoupledPressureValue.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceCoupledPressureValue_H
#define interfaceCoupledPressureValue_H

#include "patchCoupleManager.H"
#include "regionInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class interfaceCoupledPressureValue Declaration
\*---------------------------------------------------------------------------*/

class interfaceCoupledPressureValue
:
    public fixedValueFvPatchScalarField,
    public patchCoupleManager
{
    // Private data

        //- Diffusion coefficient
	    word kName_;
	    
	    //- Neighbour names
        word neighbourRegionName_;

        word neighbourPatchName_;

        //- Relaxation factor for pressure value jump
        scalar relax_;

        //- NonOrthogonal correction
        Switch nonOrthCorr_;

        //- Second order accuracy correction
        Switch secondOrder_;


    // Private functions

        //- Calculate pressure jump
        tmp<scalarField> valueJump() const;

public:

    //- Runtime type information
    TypeName("interfaceCoupledPressureValue");


    // Constructors

        //- Construct from patch and internal field
        interfaceCoupledPressureValue
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        interfaceCoupledPressureValue
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given interfaceCoupledPressureValue
        //  onto a new patch
        interfaceCoupledPressureValue
        (
            const interfaceCoupledPressureValue&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new interfaceCoupledPressureValue(*this)
            );
        }

        //- Construct as copy setting internal field reference
        interfaceCoupledPressureValue
        (
            const interfaceCoupledPressureValue&,
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
                new interfaceCoupledPressureValue(*this, iF)
            );
        }

    // Destructor

        virtual ~interfaceCoupledPressureValue();

    // Member functions
    
        //- Update the patch field coefficients
        virtual void updateCoeffs();

        //- Return patch-normal gradient
        virtual tmp<Field<scalar> > snGrad() const;

        //- Return the patch flux
        tmp<scalarField> flux() const;

        //- Return the maximum normalized coupled patch residual
        scalarField residual() const;
        
        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<scalar> > gradientBoundaryCoeffs() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
