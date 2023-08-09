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
    monolithicTemperatureFvPatchScalarField

Description
    Region couple patch field taking radiation into account

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved

SourceFiles
    monolithicTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef monolithicTemperatureFvPatchScalarField_H
#define monolithicTemperatureFvPatchScalarField_H

#include "monolithicCouplingFvPatchFields.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class monolithicTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class monolithicTemperatureFvPatchScalarField
:
    public monolithicCouplingFvPatchScalarField
{
    // Private data

        //- Name of thermal diffusivity field
        word kName_;

        //- Take radiation into account for this patch
        Switch radiation_;


public:

    //- Runtime type information
    TypeName("monolithicTemperature");


    // Constructors

        //- Construct from patch and internal field
        monolithicTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        monolithicTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given monolithicTemperatureFvPatchScalarField
        //  onto a new patch
        monolithicTemperatureFvPatchScalarField
        (
            const monolithicTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        monolithicTemperatureFvPatchScalarField
        (
            const monolithicTemperatureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new monolithicTemperatureFvPatchScalarField(*this)
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
                new monolithicTemperatureFvPatchScalarField(*this, iF)
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

            //- Return radiation
            const Switch& radiation() const
            {
                return radiation_;
            }

            //- Return access to radiation
            Switch& radiation()
            {
                return radiation_;
            }

            //- Return shadow patch field
            virtual const monolithicTemperatureFvPatchScalarField&
                shadowPatchField() const;

            //- Return neighbour field given internal cell data
            virtual tmp<scalarField> patchNeighbourField() const;

            //- Return wall temperature
            virtual tmp<scalarField> Tw() const;

            //- Return face-cell temperature
            virtual tmp<scalarField> Tc() const;


        // Evaluation functions

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

            //- Evaluate the patch field
            virtual void evaluate(const Pstream::commsTypes commsType);

            //- Evaluate the patch field
            virtual void updateCoeffs();

            //- Calculate source term
            virtual tmp<scalarField> source() const;

            //- Manipulate matrix
            virtual void manipulateMatrix(fvScalarMatrix& matrix);


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
