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

#include "regionCoupleHeatFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "regionCoupleTemperatureFvPatchScalarField.H"
#include "patchToPatchInterpolation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupleHeatFluxFvPatchScalarField::
regionCoupleHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    coupleManager_(p)
{}


Foam::regionCoupleHeatFluxFvPatchScalarField::
regionCoupleHeatFluxFvPatchScalarField
(
    const regionCoupleHeatFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    coupleManager_(ptf.coupleManager_)
{}


Foam::regionCoupleHeatFluxFvPatchScalarField::
regionCoupleHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    coupleManager_(p, dict)
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


Foam::regionCoupleHeatFluxFvPatchScalarField::
regionCoupleHeatFluxFvPatchScalarField
(
    const regionCoupleHeatFluxFvPatchScalarField& whftcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(whftcsf, iF),
    coupleManager_(whftcsf.coupleManager_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Update the patch field coefficients
void Foam::regionCoupleHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Calculate interpolated patch flux
    const polyPatch& ownPatch = patch().patch();
    const polyPatch& nbrPatch = coupleManager_.neighbourPatch().patch();
    const fvPatchField<scalar>& Tnbr = coupleManager_.neighbourPatchField<scalar>();

    patchToPatchInterpolation interpolator(nbrPatch, ownPatch);

    scalarField qNbr2own =
        interpolator.faceInterpolate
        (
            refCast<const regionCoupleTemperatureFvPatchScalarField>
            (Tnbr).flux()
        );

    // Lookup transport properties
    const dictionary& transportProperties =
        db().lookupObject<IOdictionary>("transportProperties");

    dimensionedScalar k(transportProperties.lookup("k"));

    // Enforce flux matching
    gradient() = (qNbr2own/k.value())*(-1);

    fixedGradientFvPatchScalarField::updateCoeffs();
}


//- Return the maximum normalized coupled patch residual
Foam::scalar
Foam::regionCoupleHeatFluxFvPatchScalarField::maxResidual() const
{
    // Calculate interpolated patch flux
    const polyPatch& ownPatch = patch().patch();
    const polyPatch& nbrPatch = coupleManager_.neighbourPatch().patch();
    const fvPatchField<scalar>& Tnbr = coupleManager_.neighbourPatchField<scalar>();

    patchToPatchInterpolation interpolator(nbrPatch, ownPatch);

    scalarField qNbr2own =
        interpolator.faceInterpolate
        (
            refCast<const regionCoupleTemperatureFvPatchScalarField>
            (Tnbr).flux()
        );

    // Lookup transport properties
    const dictionary& transportProperties =
        db().lookupObject<IOdictionary>("transportProperties");

    dimensionedScalar k(transportProperties.lookup("k"));

    // Calculate the maximum normalized residual
    const fvPatchScalarField& Town = *this;
    const scalarField& qOwn = k.value()*Town.snGrad();
    scalar residual =
        gMax
        (
            mag(mag(qOwn) - mag(qNbr2own))/
            max(min(gMax(mag(qOwn)),gMax(mag(qNbr2own))), SMALL)
        );

    return residual;
}


//- Write
void Foam::regionCoupleHeatFluxFvPatchScalarField::write
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
    regionCoupleHeatFluxFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
