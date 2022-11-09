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
    interfaceCoupledSpeciesValue

Author
    Constantin Habes (constantin.habes@stud.tu-darmstadt.de)

Description
    Dirichlet boundary condition for species within the ALE interface tracking.
    species jump across the interface is accounted for by adding the jump to 
    the neighbour patch value. Normal species flux is calculated and divided 
    by the density.

    Example usage:
        interfacePatchShadow
            {
                type                    interfaceCoupledSpeciesValue;
	            k			            "D";
	            neighbourRegionName     gas;
                neighbourPatchName      interfacePatch;
                nonOrthCorr             yes;
                secondOrder             yes;
                value                   uniform 0;
            }


SourceFiles
    interfaceCoupledSpeciesValue.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceCoupledSpeciesValue_H
#define interfaceCoupledSpeciesValue_H

#include "genericRegionCoupledJumpFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class interfaceCoupledSpeciesValue Declaration
\*---------------------------------------------------------------------------*/

class interfaceCoupledSpeciesValue
:
    public genericRegionCoupledJumpFvPatchField<scalar>
{

public:

    //- Runtime type information
    TypeName("interfaceCoupledSpeciesValue");


    // Constructors

        //- Construct from patch and internal field
        interfaceCoupledSpeciesValue
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        interfaceCoupledSpeciesValue
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given interfaceCoupledSpeciesValue
        //  onto a new patch
        interfaceCoupledSpeciesValue
        (
            const interfaceCoupledSpeciesValue&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new interfaceCoupledSpeciesValue(*this)
            );
        }

        //- Construct as copy setting internal field reference
        interfaceCoupledSpeciesValue
        (
            const interfaceCoupledSpeciesValue&,
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
                new interfaceCoupledSpeciesValue(*this, iF)
            );
        }

    // Member functions

        //- Calculate species jump
        virtual tmp<scalarField> valueJump() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
