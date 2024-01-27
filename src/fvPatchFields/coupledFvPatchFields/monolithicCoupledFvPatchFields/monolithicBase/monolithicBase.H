/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    monolithicBase

Description
    Region couple patch field taking radiation into account

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved

SourceFiles
    monolithicBase.C

\*---------------------------------------------------------------------------*/

#ifndef monolithicBase_H
#define monolithicBase_H

#include "monolithicCouplingFvPatchFields.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class monolithicTemperatureFvPatchScalarField;

/*---------------------------------------------------------------------------*\
                      Class monolithicBase Declaration
\*---------------------------------------------------------------------------*/

class monolithicBase
:
    public monolithicCouplingFvPatchScalarField
{
public:

    //- Runtime type information
    TypeName("monolithicBase");


    // Constructors

        //- Construct from patch and internal field
        monolithicBase
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        monolithicBase
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given monolithicBase onto a new patch
        monolithicBase
        (
            const monolithicBase&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        monolithicBase
        (
            const monolithicBase&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new monolithicBase(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new monolithicBase(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return false to disable interpolation on this patch
            virtual bool coupled() const
            {
                return false;
            }

            //- Return shadow patch field
            const monolithicBase& shadowPatchField() const;

            //- Return original patch field
            const scalarField& originalPatchField() const
            {
                return monolithicCouplingFvPatchScalarField::originalPatchField();
            }

            virtual tmp<scalarField> forig() const;

            virtual tmp<scalarField> korig() const;

            virtual tmp<scalarField> kw() const;


        // Evaluation functions

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

            //- Calculate the thermal diffusivity
            virtual tmp<scalarField> calcThermalDiffusivity
            (
                const monolithicBase& owner,
                const monolithicBase& neighbour,
                const monolithicTemperatureFvPatchScalarField& TwOwn
            ) const;

            //- Calculate the temperature
            virtual tmp<scalarField> calcTemperature
            (
                const monolithicTemperatureFvPatchScalarField& TwOwn,
                const monolithicTemperatureFvPatchScalarField& neighbour,
                const monolithicBase& ownerK
            ) const;

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix& m,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix&,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
