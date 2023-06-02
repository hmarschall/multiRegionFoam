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
    regionCoupledTemperatureJump

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
                type                    regionCoupledTemperatureJump;
	            k			            "D";
	            neighbourRegionName     gas;
                neighbourPatchName      interfacePatch;
                nonOrthCorr             yes;
                secondOrder             yes;
                value                   uniform 0;
            }


SourceFiles
    regionCoupledTemperatureJump.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledTemperatureJump_H
#define regionCoupledTemperatureJump_H

#include "genericRegionCoupledJumpFvPatchField.H"
#include "heatTransferInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class regionCoupledTemperatureJump Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledTemperatureJump
:
    public genericRegionCoupledJumpFvPatchField<scalar>
{

public:

    //- Runtime type information
    TypeName("regionCoupledTemperatureJump");


    // Constructors

        //- Construct from patch and internal field
        regionCoupledTemperatureJump
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regionCoupledTemperatureJump
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given regionCoupledTemperatureJump
        //  onto a new patch
        regionCoupledTemperatureJump
        (
            const regionCoupledTemperatureJump&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new regionCoupledTemperatureJump(*this)
            );
        }

        //- Construct as copy setting internal field reference
        regionCoupledTemperatureJump
        (
            const regionCoupledTemperatureJump&,
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
                new regionCoupledTemperatureJump(*this, iF)
            );
        }

    // Member functions

        //- Calculate species jump
        virtual tmp<scalarField> valueJump() const;

        const regionInterfaces::heatTransferInterface& heatTransInterface() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
