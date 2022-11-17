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

\*---------------------------------------------------------------------------*/

#include "regionCoupledScalarJump.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupledScalarJump::
regionCoupledScalarJump
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    genericRegionCoupledJumpFvPatchField<scalar>(p, iF)
{}


Foam::regionCoupledScalarJump::
regionCoupledScalarJump
(
    const regionCoupledScalarJump& icpv,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    genericRegionCoupledJumpFvPatchField<scalar>(icpv, p, iF, mapper)
{}


Foam::regionCoupledScalarJump::
regionCoupledScalarJump
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    genericRegionCoupledJumpFvPatchField<scalar>(p, iF, dict)
{}


Foam::regionCoupledScalarJump::
regionCoupledScalarJump
(
    const regionCoupledScalarJump& icpv,
    const DimensionedField<scalar, volMesh>& iF
)
:
    genericRegionCoupledJumpFvPatchField<scalar>(icpv, iF)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<scalarField> regionCoupledScalarJump::valueJump() const
{
    // Lookup neighbouring patch field
    const GeometricField<scalar, fvPatchField, volMesh>& nbrField = 
        nbrMesh().lookupObject<GeometricField<scalar, fvPatchField, volMesh>>
        (
            // same field name as on this side
            this->dimensionedInternalField().name()
        );

    // Calculate interpolated patch field
    Field<scalar> fieldNbrToOwn = interpolateFromNbrField<scalar>
    (
        nbrPatch()
        .patchField<GeometricField<scalar, fvPatchField, volMesh>, scalar>(nbrField)
    );

    // Get the diffusivity
    scalarField k(this->patch().size(), pTraits<scalar>::zero);

    if ( this->db().objectRegistry::foundObject<volScalarField>(kName_) )
    {
        k = this->patch().template lookupPatchField<volScalarField, scalar>(kName_);
    }
    else
    {
        k = dimensionedScalar
        (
            this->db().objectRegistry::
            lookupObject<IOdictionary>("transportProperties")
            .subDict(refPatch().boundaryMesh().mesh().name()).lookup(kName_)
        ).value();
    }

    return (fieldNbrToOwn * (k - 1));
}
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    regionCoupledScalarJump
);

} // End namespace Foam

// ************************************************************************* //
