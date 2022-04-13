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
    regionCoupleHeatFluxFvPatchScalarField

Description
    Region couple boundary condition for heat flux. Region coupling at the
    multi-region interface is accomplished via patchToPatchInterpolation,
    which is non-conservative for non-conformal boundary patches.

    Example usage:
        interfacePatch
        {
            type                    regionCoupleHeatFlux;
            neighbourRegionName     neighbourRegion;
            neighbourPatchName      neighbourPatch;
            neighbourFieldName      T;
            value                   uniform 0;
        }

SourceFiles
    regionCoupleHeatFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupleHeatFluxFvPatchScalarField_H
#define regionCoupleHeatFluxFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "coupleManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class regionCoupleHeatFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class regionCoupleHeatFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Couple manager object
        coupleManager coupleManager_;


public:

    //- Runtime type information
    TypeName("regionCoupleHeatFlux");


    // Constructors

        //- Construct from patch and internal field
        regionCoupleHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regionCoupleHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given regionCoupleHeatFluxFvPatchScalarField
        //  onto a new patch
        regionCoupleHeatFluxFvPatchScalarField
        (
            const regionCoupleHeatFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new regionCoupleHeatFluxFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        regionCoupleHeatFluxFvPatchScalarField
        (
            const regionCoupleHeatFluxFvPatchScalarField&,
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
                new regionCoupleHeatFluxFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the patch field coefficients
        virtual void updateCoeffs();

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
