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
    regionCoupleFvPatchScalarField

Description
    Region couple jump boundary condition for size Psi. Region coupling at the 
    multi-region interface is accomplished via inverse-distance weighted scalar
    interpolation and includes the possibility of a "jump" or "resistive"
    condition across the interface (e.g., T_region1 = K*T_region2).

    Example usage:
        interfacePatch
        {
            type                    regionCoupleJump;
            neighbourRegionName     neighbourRegion;
            neighbourPatchName      neighbourPatch;
            neighbourFieldName      Psi;
            k			    kName;
            value                   uniform 0;
        }

SourceFiles
    regionCoupleJumpFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupleJumpFvPatchScalarField_H
#define regionCoupleJumpFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "coupleManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class regionCoupleJumpFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class regionCoupleJumpFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Couple manager object
        coupleManager coupleManager_;

        //- "Jump" or "resistive" condition
        word kName_;


public:

    //- Runtime type information
    TypeName("regionCoupleJump");


    // Constructors

        //- Construct from patch and internal field
        regionCoupleJumpFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regionCoupleJumpFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given regionCoupleJumpFvPatchScalarField
        //  onto a new patch
        regionCoupleJumpFvPatchScalarField
        (
            const regionCoupleJumpFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new regionCoupleJumpFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        regionCoupleJumpFvPatchScalarField
        (
            const regionCoupleJumpFvPatchScalarField&,
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
                new regionCoupleJumpFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

	// Access functions
	//- Return kName
	const word& kName() const
	{
	    return kName_;
	}

	//- Return access to kName
	word& kName()
	{
	    return kName_;
	}

	// Evaluation functions
        //- Update the patch field coefficients
        virtual void updateCoeffs();

        //- Return the patch flux
        tmp<scalarField> flux() const;

        //- Return the maximum normalized coupled patch residual
        scalar maxResidual() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
