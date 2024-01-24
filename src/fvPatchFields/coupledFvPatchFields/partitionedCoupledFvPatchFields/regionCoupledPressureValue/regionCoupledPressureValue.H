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
    regionCoupledPressureValue

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
                type                    regionCoupledPressureValue;
	            k			            "v";
	            neighbourRegionName     gas;
                neighbourPatchName      interfacePatch;
                nonOrthCorr             yes;
                secondOrder             yes;
                value                   uniform 0;
            }


SourceFiles
    regionCoupledPressureValue.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledPressureValue_H
#define regionCoupledPressureValue_H

#include "genericRegionCoupledJumpFvPatchField.H"
#include "capillaryInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class regionCoupledPressureValue Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledPressureValue
:
    public genericRegionCoupledJumpFvPatchField<scalar>
{

public:

    //- Runtime type information
    TypeName("regionCoupledPressureValue");


    // Constructors

        //- Construct from patch and internal field
        regionCoupledPressureValue
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regionCoupledPressureValue
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given regionCoupledPressureValue
        //  onto a new patch
        regionCoupledPressureValue
        (
            const regionCoupledPressureValue&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new regionCoupledPressureValue(*this)
            );
        }

        //- Construct as copy setting internal field reference
        regionCoupledPressureValue
        (
            const regionCoupledPressureValue&,
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
                new regionCoupledPressureValue(*this, iF)
            );
        }

    // Member functions

        //- Calculate pressure jump
        virtual tmp<scalarField> valueJump() const;

        const regionInterfaces::capillaryInterface& capInterface() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
