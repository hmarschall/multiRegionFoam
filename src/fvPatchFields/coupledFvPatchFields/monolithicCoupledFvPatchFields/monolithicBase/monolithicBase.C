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

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "monolithicTemperatureFvPatchScalarField.H"
#include "monolithicBase.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvMatrices.H"
#include "magLongDelta.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::monolithicBase::monolithicBase
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    monolithicCouplingFvPatchScalarField(p, iF)
{}


Foam::monolithicBase::monolithicBase
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    monolithicCouplingFvPatchScalarField(p, iF, dict)
{}


Foam::monolithicBase::monolithicBase
(
    const monolithicBase& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    monolithicCouplingFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::monolithicBase::monolithicBase
(
    const monolithicBase& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    monolithicCouplingFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return a named shadow patch field
const Foam::monolithicBase&
Foam::monolithicBase::shadowPatchField() const
{
    return dynamic_cast<const monolithicBase&>
    (
        monolithicCouplingFvPatchScalarField::shadowPatchField()
    );
}


// Return a named shadow patch field
Foam::tmp<Foam::scalarField> Foam::monolithicBase::forig() const
{
    if
    (
        dimensionedInternalField().dimensions()
     == dimensionSet(1, -1, -1, 0, 0)
    )
    {
        const fvPatch& p = patch();

        const basicThermo& thermo = db().lookupObject<basicThermo>
        (
            "thermophysicalProperties"
        );

        const scalarField Tw =
            p.lookupPatchField<volScalarField, scalar>("T");

        return originalPatchField()*thermo.Cp(Tw, p.index());
    }
    else
    {
        return originalPatchField();
    }
}


// Return a named shadow patch field
Foam::tmp<Foam::scalarField> Foam::monolithicBase::korig() const
{
    const fvPatch& p = patch();
    const magLongDelta& mld = magLongDelta::New(p.boundaryMesh().mesh());

    return forig()/(1 - p.weights())/mld.magDelta(p.index());
}


// Return a named shadow patch field
Foam::tmp<Foam::scalarField> Foam::monolithicBase::kw() const
{
    if
    (
        dimensionedInternalField().dimensions()
     == dimensionSet(1, -1, -1, 0, 0, 0, 0)
    )
    {
        const fvPatch& p = patch();

        const basicThermo& thermo = db().lookupObject<basicThermo>
        (
            "thermophysicalProperties"
        );

        const scalarField Tw =
            p.lookupPatchField<volScalarField, scalar>("T");

        return *this*thermo.Cp(Tw, p.index());
    }
    else
    {
        return *this;
    }
}


Foam::tmp<Foam::scalarField>
Foam::monolithicBase::calcThermalDiffusivity
(
    const monolithicBase& owner,
    const monolithicBase& neighbour,
    const monolithicTemperatureFvPatchScalarField& TwOwn
) const
{
    return shadowPatchField().calcThermalDiffusivity(owner, neighbour, TwOwn);
}


Foam::tmp<Foam::scalarField>
Foam::monolithicBase::calcTemperature
(
    const monolithicTemperatureFvPatchScalarField& TwOwn,
    const monolithicTemperatureFvPatchScalarField& neighbour,
    const monolithicBase& ownerK
) const
{
    return shadowPatchField().calcTemperature(TwOwn, neighbour, ownerK);
}


void Foam::monolithicBase::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    updateCoeffs();
}


void
Foam::monolithicBase::initInterfaceMatrixUpdate
(
    const scalarField&,
    scalarField&,
    const lduMatrix&,
    const scalarField&,
    const direction,
    const Pstream::commsTypes,
    const bool switchToLhs
) const
{
    FatalErrorIn
    (
        "monolithicBase::initInterfaceMatrixUpdate"
    )   << "Undefined function: this patch field cannot be used "
        << "on active variables"
        << abort(FatalError);
}


void
Foam::monolithicBase::updateInterfaceMatrix
(
    const scalarField&,
    scalarField&,
    const lduMatrix&,
    const scalarField&,
    const direction,
    const Pstream::commsTypes,
    const bool switchToLhs
) const
{
    FatalErrorIn
    (
        "monolithicBase::updateInterfaceMatrix"
    )   << "Undefined function: this patch field cannot be used "
        << "on active variables"
        << abort(FatalError);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    monolithicBase
);

} // End namespace Foam

// ************************************************************************* //
