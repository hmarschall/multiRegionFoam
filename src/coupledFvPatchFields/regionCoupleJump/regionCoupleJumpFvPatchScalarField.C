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

#include "regionCoupleJumpFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "patchToPatchInterpolation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupleJumpFvPatchScalarField::
regionCoupleJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    coupleManager_(p),
    kName_("none")
{}


Foam::regionCoupleJumpFvPatchScalarField::
regionCoupleJumpFvPatchScalarField
(
    const regionCoupleJumpFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    coupleManager_(ptf.coupleManager_),
    kName_(ptf.kName_)
{}


Foam::regionCoupleJumpFvPatchScalarField::
regionCoupleJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
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


Foam::regionCoupleJumpFvPatchScalarField::
regionCoupleJumpFvPatchScalarField
(
    const regionCoupleJumpFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wtcsf, iF),
    coupleManager_(wtcsf.coupleManager_),
    kName_(wtcsf.kName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Update the patch field coefficients
void Foam::regionCoupleJumpFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Calculate interpolated patch field
    const polyPatch& ownPatch = patch().patch();
    const polyPatch& nbrPatch = coupleManager_.neighbourPatch().patch();
    const fvPatchField<scalar>& psiNbr = coupleManager_.neighbourPatchField<scalar>();

    patchToPatchInterpolation interpolator(nbrPatch, ownPatch);
    scalarField psiNbr2Own = interpolator.faceInterpolate(psiNbr);

    // Lookup coupled solution controls
//    const dictionary& coupledSolutionDict =
//        db().time().lookupObject<IOdictionary>("coupledSolutionDict");

//    const scalar& relax =
//        readScalar(coupledSolutionDict
//        .subDict("partitioned").lookup("coupleRelaxFactor"));

    scalar relax = 0.2;

    // Lookup diffusivity field
    const fvPatchScalarField& kpf =
	lookupPatchField<volScalarField, scalar>(kName_);
    // Enforce psi jump boundary condition
    operator==(*this + relax*(kpf*psiNbr2Own - *this));

    fixedValueFvPatchScalarField::updateCoeffs();
}


//- Return the patch flux
Foam::tmp<Foam::scalarField>
Foam::regionCoupleJumpFvPatchScalarField::flux() const
{
    const fvPatchScalarField& kpf =
	lookupPatchField<volScalarField, scalar>(kName_);

    const fvPatchScalarField& psi = *this;

    return kpf*psi.snGrad();
}


//- Return the maximum normalized coupled patch residual
Foam::scalar
Foam::regionCoupleJumpFvPatchScalarField::maxResidual() const
{
    // Calculate interpolated patch field
    const polyPatch& ownPatch = patch().patch();
    const polyPatch& nbrPatch = coupleManager_.neighbourPatch().patch();
    const fvPatchField<scalar>& psiNbr = coupleManager_.neighbourPatchField<scalar>();

    patchToPatchInterpolation interpolator(nbrPatch, ownPatch);
    scalarField psiNbr2Own = interpolator.faceInterpolate(psiNbr);

    // Lookup diffusivity field
    const fvPatchScalarField& kpf =
	lookupPatchField<volScalarField, scalar>(kName_);

    // Calculate the maximum normalized residual
    const scalarField& psiOwn = *this;
    scalar residual =
        gMax
        (
            mag(psiOwn - kpf*psiNbr2Own)/
            max(min(gMax(psiOwn),gMax(kpf*psiNbr2Own)), SMALL)
        );

    return residual;
}


//- Write
void Foam::regionCoupleJumpFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    coupleManager_.writeEntries(os);
    writeEntry("k", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    regionCoupleJumpFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
