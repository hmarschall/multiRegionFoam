/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    twoPhaseInterfaceVelocity

Description
    Generalized grid interface patch field.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

SourceFiles
    twoPhaseInterfaceVelocity.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseInterfaceVelocity_H
#define twoPhaseInterfaceVelocity_H

#include "ggiFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseInterfaceVelocity Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseInterfaceVelocityFvPatchVectorField
:
    public ggiFvPatchField<vector>
{

    // Private data

        //- Whether this patch is master
        Switch master_;

        //- Viscous force jump
        vectorField forceJump_;

    // Private functions

public:

    //- Runtime type information
    TypeName("twoPhaseInterfaceVelocity");


    // Constructors

        //- Construct from patch and internal field
        twoPhaseInterfaceVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        twoPhaseInterfaceVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given twoPhaseInterfaceVelocity 
        //  onto a new patch
        twoPhaseInterfaceVelocityFvPatchVectorField
        (
            const twoPhaseInterfaceVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        twoPhaseInterfaceVelocityFvPatchVectorField
        (
            const twoPhaseInterfaceVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new twoPhaseInterfaceVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new twoPhaseInterfaceVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return true if master patch
            const Switch& master() const;

            //- Return reference to forceJump at the interface
            const vectorField& forceJump() const
            {
                return forceJump_;
            }

        // Evaluation functions

//             //- Return neighbour field given internal cell data
//             virtual tmp<Field<vector> > patchNeighbourField() const;

//             //- Update the coefficients associated with the patch field
//             virtual void updateCoeffs();

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

//             //- Evaluate the patch field
//             virtual void evaluate
//             (
//                 const Pstream::commsTypes commsType
//             );

//             //- Transform neighbour field
//             virtual void transformCoupleField
//             (
//                 vectorField& f,
//                 const direction cmpt
//             ) const
//             {
//                 ggiLduInterfaceField::transformCoupleField(f, cmpt);
//             }

//             //- Initialise neighbour matrix update
//             virtual void initInterfaceMatrixUpdate
//             (
//                 const scalarField& psiInternal,
//                 scalarField& result,
//                 const lduMatrix& m,
//                 const scalarField& coeffs,
//                 const direction cmpt,
//                 const Pstream::commsTypes commsType,
//                 const bool switchToLhs
//             ) const;

//             //- Update result field based on interface functionality
//             virtual void updateInterfaceMatrix
//             (
//                 const scalarField& psiInternal,
//                 scalarField& result,
//                 const lduMatrix&,
//                 const scalarField& coeffs,
//                 const direction cmpt,
//                 const Pstream::commsTypes commsType,
//                 const bool switchToLhs
//             ) const;

            //- Calculate patch face interpolate given weighting factors
            virtual void patchInterpolate
            (
                GeometricField<vector, fvsPatchField, surfaceMesh>& fField,
                const scalarField& pL
            ) const;

            //- Calculate patch face interpolate given two weighting factors
            virtual void patchInterpolate
            (
                GeometricField<vector, fvsPatchField, surfaceMesh>& fField,
                const scalarField& pL,
                const scalarField& pY
            ) const;

            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<vector>& matrix);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
