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

#include "regionCoupleTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "patchToPatchInterpolation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupleTemperatureFvPatchScalarField::
regionCoupleTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    coupleManager_(p),
    K_(0)
{}


Foam::regionCoupleTemperatureFvPatchScalarField::
regionCoupleTemperatureFvPatchScalarField
(
    const regionCoupleTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    coupleManager_(ptf.coupleManager_),
    K_(ptf.K_)
{}


Foam::regionCoupleTemperatureFvPatchScalarField::
regionCoupleTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    coupleManager_(p, dict),
    K_(readScalar(dict.lookup("K")))
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


Foam::regionCoupleTemperatureFvPatchScalarField::
regionCoupleTemperatureFvPatchScalarField
(
    const regionCoupleTemperatureFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wtcsf, iF),
    coupleManager_(wtcsf.coupleManager_),
    K_(wtcsf.K_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Update the patch field coefficients
void Foam::regionCoupleTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Calculate interpolated patch field
    const polyPatch& ownPatch = patch().patch();
    const polyPatch& nbrPatch = coupleManager_.neighbourPatch().patch();
    const fvPatchField<scalar>& Tnbr = coupleManager_.neighbourPatchField<scalar>();

    patchToPatchInterpolation interpolator(nbrPatch, ownPatch);
    scalarField Tnbr2own = interpolator.faceInterpolate(Tnbr);

    // Lookup coupled solution controls
//    const dictionary& coupledSolutionDict =
//        db().time().lookupObject<IOdictionary>("coupledSolutionDict");

//    const scalar& relax =
//        readScalar(coupledSolutionDict
//        .subDict("partitioned").lookup("coupleRelaxFactor"));

    scalar relax = 0.5;

    // Enforce temperature boundary condition
    operator==(*this + relax*(K_*Tnbr2own - *this));

    fixedValueFvPatchScalarField::updateCoeffs();
}


//- Return the patch flux
Foam::tmp<Foam::scalarField>
Foam::regionCoupleTemperatureFvPatchScalarField::flux() const
{
    const dictionary& transportProperties =
        db().lookupObject<IOdictionary>("transportProperties");

    dimensionedScalar k(transportProperties.lookup("k"));

    const fvPatchScalarField& T = *this;

    return k.value()*T.snGrad();
}


//- Return the maximum normalized coupled patch residual
Foam::scalar
Foam::regionCoupleTemperatureFvPatchScalarField::maxResidual() const
{
    // Calculate interpolated patch field
    const polyPatch& ownPatch = patch().patch();
    const polyPatch& nbrPatch = coupleManager_.neighbourPatch().patch();
    const fvPatchField<scalar>& Tnbr = coupleManager_.neighbourPatchField<scalar>();

    patchToPatchInterpolation interpolator(nbrPatch, ownPatch);
    scalarField Tnbr2own = interpolator.faceInterpolate(Tnbr);

    // Calculate the maximum normalized residual
    const scalarField& Town = *this;
    scalar residual =
        gMax
        (
            mag(Town - K_*Tnbr2own)/
            max(min(gMax(Town),gMax(K_*Tnbr2own)), SMALL)
        );

    return residual;
}


//- Write
void Foam::regionCoupleTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    coupleManager_.writeEntries(os);
    os.writeKeyword("K") << K_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    regionCoupleTemperatureFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //