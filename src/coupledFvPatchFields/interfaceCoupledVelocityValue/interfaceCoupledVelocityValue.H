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
    interfaceCoupledVelocityValue

Author
    Holger Marschall (marschall@csi.tu-darmstadt.de)

Description
    Couple boundary condition for velocity within the
    ALE interface tracking framework.
    The no-jump condition of the velocity at the interface is taken into
    account by copying the value from the neighbour patch.
    The normal flux is calculated as the normal gradient multiplied by mu.

    Example usage:
        interfacePatch
        {
            type                    interfaceCoupledVelocityValue;
            neighbourPatchName      neighbourPatch;
            k                       "mu";
            value                   uniform (0,0,0);
        }

SourceFiles
    interfaceCoupledVelocityValue.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceCoupledVelocityValue_H
#define interfaceCoupledVelocityValue_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "fixedValueFvPatchFields.H"

#include "patchCoupleManager.H"
//#include "interfacialTransport.H"
#include "surfaceTracking.H"
#include "meshInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class interfaceCoupledVelocityValue Declaration
\*---------------------------------------------------------------------------*/

class interfaceCoupledVelocityValue
:
    public fixedValueFvPatchVectorField //fixedValueFvPatchVectorField
{
    // Private data

        //- Diffusion coefficient
	    word k_;

        //relaxation factor for velocity jump
        scalar relax_;

        //- Name of flux field
        word phiName_;

        //- Name of density field
        word rhoName_;

        //- NonOrthogonal correction
        Switch nonOrthCorr_;

        //- Second order accuracy correction
        Switch secondOrder_;

        //- Current time index
        label curTimeIndex_;

        // Old face centres of patch
        // for (second order) time integration
        vectorField Fc_;
        vectorField oldFc_;
        vectorField oldOldFc_;

        //- Couple manager
        autoPtr<patchCoupleManager> coupleManagerPtr_;

        //- Interfacial transport
        mutable autoPtr<interfacialTransport> interfacePtr_;

    // Private functions

        //- Reference to meshInterface
        const meshInterface& ale() const;

        //- Calc pressure jump
        tmp<vectorField> curJump() const;

public:

    //- Runtime type information
    TypeName("interfaceCoupledVelocityValue");


    // Constructors

        //- Construct from patch and internal field
        interfaceCoupledVelocityValue
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        interfaceCoupledVelocityValue
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given interfaceCoupledVelocityValue
        //  onto a new patch
        interfaceCoupledVelocityValue
        (
            const interfaceCoupledVelocityValue&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new interfaceCoupledVelocityValue(*this)
            );
        }

        //- Construct as copy setting internal field reference
        interfaceCoupledVelocityValue
        (
            const interfaceCoupledVelocityValue&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new interfaceCoupledVelocityValue(*this, iF)
            );
        }


    // Member functions
        //- Update the patch field coefficients
        virtual void updateCoeffs();

        //- Return patch-normal gradient
        virtual tmp<Field<vector> > snGrad() const;

        //- Return the patch flux
        tmp<vectorField> flux() const;

        //- Return the maximum normalized coupled patch residual
        scalarField residual() const;

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the gradient of this patchField
//        virtual tmp<Field<vector> > gradientInternalCoeffs() const;
        
        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the gradient of this patchField
//        virtual tmp<Field<vector> > gradientBoundaryCoeffs() const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<vector> > gradientBoundaryCoeffs() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

