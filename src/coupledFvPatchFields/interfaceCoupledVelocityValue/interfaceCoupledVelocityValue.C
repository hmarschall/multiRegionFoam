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

#include "interfaceCoupledVelocityValue.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "fvCFD.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCoupledVelocityValue::
interfaceCoupledVelocityValue
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    patchCoupleManager(p),
    kName_("k"),
    relax_(1.0),
    neighbourRegionName_(),
    neighbourPatchName_(),
    phiName_("phi"),
    rhoName_("rhoA"),
    nonOrthCorr_(false),
    secondOrder_(false)
{}


Foam::interfaceCoupledVelocityValue::
interfaceCoupledVelocityValue
(
    const interfaceCoupledVelocityValue& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    patchCoupleManager(p),
    kName_(ptf.kName_),
    relax_(ptf.relax_),
    neighbourRegionName_(ptf.neighbourRegionName_),
    neighbourPatchName_(ptf.neighbourPatchName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    nonOrthCorr_(ptf.nonOrthCorr_),
    secondOrder_(ptf.secondOrder_)
{}


Foam::interfaceCoupledVelocityValue::
interfaceCoupledVelocityValue
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    patchCoupleManager(p, dict),
    kName_(dict.lookupOrDefault<word>("k", word::null)),
    relax_(dict.lookupOrDefault<scalar>("relax",1.0)),
    neighbourRegionName_
    (
        dict.lookupOrDefault<word>("neighbourRegionName", word::null)
    ),
    neighbourPatchName_
    (
        dict.lookupOrDefault<word>("neighbourPatchName", word::null)
    ),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    nonOrthCorr_(false),
    secondOrder_(false)
{
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        evaluate();
    }
}


Foam::interfaceCoupledVelocityValue::
interfaceCoupledVelocityValue
(
    const interfaceCoupledVelocityValue& icvv,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(icvv, iF),
    patchCoupleManager(icvv),
    kName_(icvv.kName_),
    relax_(icvv.relax_),
    neighbourRegionName_(icvv.neighbourRegionName_),
    neighbourPatchName_(icvv.neighbourPatchName_),
    phiName_(icvv.phiName_),
    rhoName_(icvv.rhoName_),
    nonOrthCorr_(icvv.nonOrthCorr_),
    secondOrder_(icvv.secondOrder_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Update the patch field coefficients
void Foam::interfaceCoupledVelocityValue::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Calculate interpolated patch field
    vectorField fieldNbrToOwn(patch().size(), vector::zero);

    // Lookup neighbouring patch field
    const volVectorField& nbrField = nbrMesh().lookupObject<volVectorField>
        (
            //presume same field name as on this side
            this->dimensionedInternalField().name()
        );

    fieldNbrToOwn = interpolateFromNbrField<vector>
        (
            nbrPatch().patchField<volVectorField, vector>(nbrField)
        );

    // Enforce flux continuity at interface

    //- Non-const access to flux on patch
    fvsPatchField<scalar>& patchPhiField = const_cast<fvsPatchField<scalar>& >
    (
        this->db().lookupObject<surfaceScalarField>(phiName_)
        .boundaryField()[this->patch().index()]
    );

    //- Get the flux on neighboring patch 
    surfaceScalarField nbrPhi =
        nbrMesh().lookupObject<surfaceScalarField>(phiName_);

    //- Impose interpolated flux field
    patchPhiField = interpolateFromNbrField<scalar>
        (
            nbrPatch().patchField<surfaceScalarField, scalar>(nbrPhi)
        )*(-1.); // consider outer normals pointing in opposite directions

    // Enforce fixed value condition
    operator==
    (
        *this 
      + relax_*
        (
            fieldNbrToOwn 
          - *this
        )
    );

    fixedValueFvPatchVectorField::updateCoeffs();
}

tmp<Field<vector> > interfaceCoupledVelocityValue::gradientBoundaryCoeffs() 
const
{
    bool secondOrder_ = false;

    word UName = this->dimensionedInternalField().name();

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + UName + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = delta - n*(n&delta);

    vectorField dUP = (k&gradU.patchInternalField());

    if (secondOrder_)
    {
        vectorField nGradUP = (n&gradU.patchInternalField());

        vectorField nGradU =
            2
           *(
                *this
              - (patchInternalField() + dUP)
            )*this->patch().deltaCoeffs()
          - nGradUP;

        vectorField nGradUn = n*(n&nGradU);

        return
            this->patch().deltaCoeffs()
           *(
                2*(*this - dUP)
              - patchInternalField()
            )
          - nGradUP
          - nGradUn;
    }

    // First order
    vectorField nGradU =
        (
            *this
          - (patchInternalField() + dUP)
        )*this->patch().deltaCoeffs();

    vectorField nGradUn = n*(n&nGradU);

    return
        this->patch().deltaCoeffs()
       *(
           *this - dUP
        )
      - nGradUn;
}

//- Return the patch-normal gradient
Foam::tmp<Foam::Field<vector> >
interfaceCoupledVelocityValue::snGrad() const
{
//    bool secondOrder_ = true;

    word UName = this->dimensionedInternalField().name();

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + UName + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = delta - n*(n&delta);

    vectorField dUP(this->patch().size(), vector::zero);

    if (nonOrthCorr_)
    {
        dUP = (k&gradU.patchInternalField());
    }

    if (secondOrder_)
    {
        vectorField nGradUP = (n&gradU.patchInternalField());

        tmp<Field<vector> > tnGradU
        (
            new vectorField(this->patch().size(), vector::zero)
        );

        tnGradU() =
            2
           *(
                *this
              - (patchInternalField() + dUP)
            )*this->patch().deltaCoeffs()
          - nGradUP;

        tnGradU() -= n*(n&tnGradU());

        return tnGradU;
    }


    // First order
    tmp<Field<vector> > tnGradU
    (
        new vectorField(this->patch().size(), vector::zero)
    );

    tnGradU() =
        (
            *this
          - (patchInternalField() + dUP)
        )*this->patch().deltaCoeffs();

    tnGradU() -= n*(n&tnGradU());

    return tnGradU;
}

//- Return the patch flux
Foam::tmp<Foam::vectorField>
Foam::interfaceCoupledVelocityValue::flux() const
{
    dimensionedScalar k
    (
        nbrMesh().lookupObject<IOdictionary>
        ("surfaceProperties").lookup(kName_)
    );

    return (this->snGrad()*k.value());
}


//- Return the maximum normalized coupled patch residual
Foam::scalarField Foam::interfaceCoupledVelocityValue::residual() const
{
    // Calculate interpolated patch field
    vectorField fieldNbrToOwn(patch().size(), vector::zero);

    // Lookup neighbouring patch field
    const volVectorField& nbrField = nbrMesh().lookupObject<volVectorField>
        (
            //presume same field name as on this side
            this->dimensionedInternalField().name()
        );

    // Interpolate veclocity face values
    fieldNbrToOwn = interpolateFromNbrField<vector>
        (
            nbrPatch().patchField<volVectorField, vector>(nbrField)
        );

    // Get the fluxes
    const fvsPatchField<scalar>& phip =
    (
        this->db().lookupObject<surfaceScalarField>(phiName_)
        .boundaryField()[this->patch().index()]
    );

    const surfaceScalarField& nbrPhi =
        nbrMesh().lookupObject<surfaceScalarField>(phiName_);

    scalarField nbrPhip = interpolateFromNbrField<scalar>
        (
            nbrPatch().patchField<surfaceScalarField, scalar>(nbrPhi)
        );

	// Calculate the maximal normalized residual
    const vectorField& fown = *this;

    const scalarField& residual =
        (
            mag(phip - nbrPhip)
           /min
            (
                mag(fown&patch().nf()),
                mag(fieldNbrToOwn&patch().nf())
            )
        );

    return residual;
}

//- Write
void Foam::interfaceCoupledVelocityValue::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("k") << kName_ << token::END_STATEMENT << nl;
    os.writeKeyword("relax") << relax_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourRegionName") << neighbourRegionName_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourPatchName") << neighbourPatchName_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("phiName") << phiName_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("rhoName") << rhoName_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("nonOrthCorr") << nonOrthCorr_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("secondOrder") << secondOrder_ 
        << token::END_STATEMENT << nl;
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

makePatchTypeField
(
    fvPatchVectorField,
    interfaceCoupledVelocityValue
);

} // End namespace Foam

// ************************************************************************* //

