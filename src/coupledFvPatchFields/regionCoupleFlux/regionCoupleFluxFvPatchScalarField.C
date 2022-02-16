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

#include "regionCoupleFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "regionCoupleJumpFvPatchScalarField.H"
#include "patchToPatchInterpolation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupleFluxFvPatchScalarField::
regionCoupleFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    coupleManager_(p),
    kName_("none")
{}


Foam::regionCoupleFluxFvPatchScalarField::
regionCoupleFluxFvPatchScalarField
(
    const regionCoupleFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    coupleManager_(ptf.coupleManager_),
    kName_(ptf.kName_)
{}


Foam::regionCoupleFluxFvPatchScalarField::
regionCoupleFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    coupleManager_(p, dict),
    kName_(dict.lookup("k"))
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        evaluate();
    }
}


Foam::regionCoupleFluxFvPatchScalarField::
regionCoupleFluxFvPatchScalarField
(
    const regionCoupleFluxFvPatchScalarField& whftcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(whftcsf, iF),
    coupleManager_(whftcsf.coupleManager_),
    kName_(whftcsf.kName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Update the patch field coefficients
void Foam::regionCoupleFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Calculate interpolated patch flux
    const polyPatch& ownPatch = patch().patch();
    const polyPatch& nbrPatch = coupleManager_.neighbourPatch().patch();
    const fvPatchField<scalar>& psiNbr = coupleManager_.neighbourPatchField<scalar>();

    patchToPatchInterpolation interpolator(nbrPatch, ownPatch);

    scalarField jNbr2Own =
        interpolator.faceInterpolate
        (
            refCast<const regionCoupleJumpFvPatchScalarField>
            (psiNbr).flux()
        );

    // Lookup diffusivity field
    const fvPatchScalarField& k =
	lookupPatchField<volScalarField, scalar>(kName_);

    // Enforce flux matching
    gradient() = (jNbr2Own/k)*(-1);

    fixedGradientFvPatchScalarField::updateCoeffs();
}


//- Return the maximum normalized coupled patch residual
Foam::scalar
Foam::regionCoupleFluxFvPatchScalarField::maxResidual() const
{
    // Calculate interpolated patch flux
    const polyPatch& ownPatch = patch().patch();
    const polyPatch& nbrPatch = coupleManager_.neighbourPatch().patch();
    const fvPatchField<scalar>& psiNbr = coupleManager_.neighbourPatchField<scalar>();

    patchToPatchInterpolation interpolator(nbrPatch, ownPatch);

    scalarField jNbr2Own =
        interpolator.faceInterpolate
        (
            refCast<const regionCoupleJumpFvPatchScalarField>
            (psiNbr).flux()
        );

    // Lookup diffusivity field
    const fvPatchScalarField& k =
	lookupPatchField<volScalarField, scalar>(kName_);

    // Calculate the maximum normalized residual
    const fvPatchScalarField& psiOwn = *this;
    const scalarField& jOwn = k*psiOwn.snGrad();
    scalar residual =
        gMax
        (
            mag(mag(jOwn) - mag(jNbr2Own))/
            max(min(gMax(mag(jOwn)),gMax(mag(jNbr2Own))), SMALL)
        );

    return residual;
}


//- Write
void Foam::regionCoupleFluxFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    coupleManager_.writeEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    regionCoupleFluxFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
