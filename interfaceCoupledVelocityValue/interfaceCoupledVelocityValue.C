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
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "fvCFD.H"
#include "ggiInterpolation.H"

#include "backwardDdtScheme.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::meshInterface& 
Foam::interfaceCoupledVelocityValue::ale() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const objectRegistry& obr = mesh.objectRegistry::parent();

    if (!obr.foundObject<meshInterface>("aleProperties"))
    {
        FatalErrorIn("interfaceCoupledVelocityValue::")
            << "meshInterface object not found but required."
            << abort(FatalError);
    }

    return obr.lookupObject<meshInterface>("aleProperties");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCoupledVelocityValue::
interfaceCoupledVelocityValue
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    k_("k"),
    relax_(1.0),
    phiName_("phi"),
    rhoName_("rhoA"),
    nonOrthCorr_(false),
    secondOrder_(false),
    curTimeIndex_( dimensionedInternalField().mesh().time().timeIndex() ),
    Fc_( p.patch().size(), vector::zero ),
    oldFc_( p.patch().size(), vector::zero ),
    oldOldFc_( p.patch().size(), vector::zero )
{
    coupleManagerPtr_.reset
    (
        new patchCoupleManager(p)
    );

    interfacePtr_.reset
    (
        new interfacialTransport
        (
//            word::null,
            this->db().lookupObject<volVectorField>("U"),
            this->db().lookupObject<surfaceScalarField>("phi"),
            this->patch().name(),
            coupleManagerPtr_->neighbourPatch().patch().boundaryMesh().mesh()
            .lookupObject<IOdictionary>("surfaceProperties")
        )
    );
}


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
    k_(ptf.k_),
    relax_(ptf.relax_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    nonOrthCorr_(ptf.nonOrthCorr_),
    secondOrder_(ptf.secondOrder_),
    curTimeIndex_(ptf.curTimeIndex_),
    Fc_(p.patch().size(), vector::zero),
    oldFc_(p.patch().size(), vector::zero),
    oldOldFc_(p.patch().size(), vector::zero)
{
    coupleManagerPtr_.reset
    (
        new patchCoupleManager
        (
            ptf.coupleManagerPtr_
        )
    );

    interfacePtr_.reset
    (
        new interfacialTransport
        (
//            word::null,
            this->db().lookupObject<volVectorField>("U"),
            this->db().lookupObject<surfaceScalarField>("phi"),
            this->patch().name(),
            coupleManagerPtr_->neighbourPatch().patch().boundaryMesh().mesh()
            .lookupObject<IOdictionary>("surfaceProperties")
        )
    );
}


Foam::interfaceCoupledVelocityValue::
interfaceCoupledVelocityValue
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    k_(dict.lookupOrDefault<word>("k", word::null)),
    relax_(dict.lookupOrDefault<scalar>("relax",1.0)),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    nonOrthCorr_(false),
    secondOrder_(false),
    curTimeIndex_(dimensionedInternalField().mesh().time().timeIndex()),
    Fc_(p.patch().size(), vector::zero),
    oldFc_(p.patch().size(), vector::zero),
    oldOldFc_(p.patch().size(), vector::zero)
{
    coupleManagerPtr_.reset
    (
        new patchCoupleManager
        (
            p, dict
        )
    );

    interfacePtr_.reset
    (
        new interfacialTransport
        (
//            word::null,
            this->db().lookupObject<volVectorField>("U"),
            this->db().lookupObject<surfaceScalarField>("phi"),
            this->patch().name(),
            coupleManagerPtr_->neighbourPatch().patch().boundaryMesh().mesh()
            .lookupObject<IOdictionary>("surfaceProperties")
        )
    );

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
    const interfaceCoupledVelocityValue& wtcsf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(wtcsf, iF),
    k_(wtcsf.k_),
    relax_(wtcsf.relax_),
    phiName_(wtcsf.phiName_),
    rhoName_(wtcsf.rhoName_),
    nonOrthCorr_(wtcsf.nonOrthCorr_),
    secondOrder_(wtcsf.secondOrder_),
    curTimeIndex_(wtcsf.curTimeIndex_),
    Fc_(wtcsf.oldFc_),
    oldFc_(wtcsf.oldFc_),
    oldOldFc_(wtcsf.oldOldFc_)
{
    coupleManagerPtr_.reset
    (
        new patchCoupleManager
        (
            wtcsf.coupleManagerPtr_
        )
    );

    interfacePtr_.reset
    (
        new interfacialTransport
        (
//            word::null,
            this->db().lookupObject<volVectorField>("U"),
            this->db().lookupObject<surfaceScalarField>("phi"),
            this->patch().name(),
            coupleManagerPtr_->neighbourPatch().patch().boundaryMesh().mesh()
            .lookupObject<IOdictionary>("surfaceProperties")
        )
    );
}


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

//    const fvMesh& mesh = patch().boundaryMesh().mesh();
//    const objectRegistry& obr = mesh.objectRegistry::parent();

//    if (!obr.foundObject<meshInterface>("aleProperties"))
//    {
//        FatalErrorIn("interfaceCoupledPressureValue::updateCoeffs(...)")
//            << "meshInterface object not found but required."
//            << abort(FatalError);
//    }

//    const meshInterface& ale =
//        obr.lookupObject<meshInterface>("aleProperties");

    // Lookup neighbouring patch field
    const volVectorField& nbrField = ale().meshA().lookupObject<volVectorField>
        (
            //presume same field name as on this side
            this->dimensionedInternalField().name()
        );

    const fvPatchVectorField& nbrPatchField =
    (
        ale().patch().patchField<volVectorField, vector>(nbrField)
    );

    // Interpolate veclocity face values A-to-B
    ale().transferFacesAToB(nbrPatchField, fieldNbrToOwn);

    // Enforce flux continuity at interface

    //- Non-const access to flux on patch
    fvsPatchField<scalar>& patchPhiField = const_cast<fvsPatchField<scalar>& >
    (
        this->db().lookupObject<surfaceScalarField>(phiName_)
        .boundaryField()[this->patch().index()]
    );

    //- Get the flux on neighboring patch 
    surfaceScalarField nbrPhi =
        ale().meshA().lookupObject<surfaceScalarField>(phiName_);

    scalarField nbrPatchPhiField =
    (
        ale().patch().patchField<surfaceScalarField, scalar>(nbrPhi)
    )*(-1.); // consider outer normals pointing in opposite directions

    //- Interpolate flux face values A-to-B
    ale().transferFacesAToB(nbrPatchPhiField, patchPhiField);

    // Enforce fixed value condition
    const surfaceScalarField& phi =
        this->db().lookupObject<surfaceScalarField>(phiName_);

//    tmp<vectorField> n = p.nf();
//    const vectorField& Sf = p.Sf();

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        operator==
        (
            *this 
          + relax_*
            (
                fieldNbrToOwn 
//              + Sf*phip
//              - n()*(n()&Up)
//              + interfacePtr_->massFluxCorr()
//                *(
//                    1./interfacePtr_->rhoB()
//                  - 1./interfacePtr_->rhoA()
//                ).value()*patch().nf()
              - *this
            )
        );
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        if (!this->db().objectRegistry::found(rhoName_))
        {
            // Rho not available, do not update
            return;
        }
//        const fvPatchField<scalar>& rhop =
//            lookupPatchField<volScalarField, scalar>(rhoName_);

        operator==
        (
            *this 
          + relax_*
            (
                fieldNbrToOwn 
//              + Sf*phip/rhop
//              - n()*(n()&Up)
//              + interfacePtr_->massFluxCorr()
//                *(
//                    1./interfacePtr_->rhoB()
//                  - 1./interfacePtr_->rhoA()
//                ).value()*patch().nf()
              - *this
            )
        );
    }
    else
    {
        FatalErrorIn
        (
            "interfaceCoupledVelocityValue::updateCoeffs()"
        )
            << "dimensions of phi are incorrect\n"
            << "    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

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
    const fvMesh& nbrMesh = 
        coupleManagerPtr_->neighbourPatch().boundaryMesh().mesh();

    dimensionedScalar k
    (
      nbrMesh.lookupObject<IOdictionary>
      ("surfaceProperties").lookup(k_)
    );

    return (this->snGrad()*k.value());
}


//- Return the maximum normalized coupled patch residual
Foam::scalarField Foam::interfaceCoupledVelocityValue::residual() const
{
    // Calculate interpolated patch field
    const polyPatch& ownPatch = patch().patch();
    const polyPatch& nbrPatch = coupleManagerPtr_->neighbourPatch().patch();
    const fvPatch& nbrFvPatch = coupleManagerPtr_->neighbourPatch();

    const fvMesh& nbrMesh = 
        coupleManagerPtr_->neighbourPatch().boundaryMesh().mesh();

    const GGIInterpolation<polyPatch, polyPatch>& interpolator = 
        coupleManagerPtr_->ggiInterpolator(nbrPatch,ownPatch);

    fvPatchVectorField nbrPatchField =
         coupleManagerPtr_->neighbourPatchField<vector>();

    vectorField fieldNbrToOwn = interpolator.masterToSlave
    (
        nbrPatchField
    );

    // Get the fluxes
    const fvsPatchField<scalar>& phip =
    (
        this->db().lookupObject<surfaceScalarField>(phiName_)
        .boundaryField()[this->patch().index()]
    );

    const surfaceScalarField& nbrPhi =
        nbrMesh.lookupObject<surfaceScalarField>(phiName_);//

    scalarField nbrPhip =
        interpolator.masterToSlave
        (
            nbrFvPatch.patchField<surfaceScalarField, scalar>(nbrPhi)
        );

	// Calculate the maximal normalized residual
    const vectorField& fown = *this;

    tmp<vectorField> n = patch().nf();
//    tmp<scalarField> magS = patch().magSf();

    const scalarField& residual =
        (
            mag(phip - nbrPhip)/
            min(mag(fown&n()),mag(fieldNbrToOwn&n()))
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
    coupleManagerPtr_->writeEntries(os);
    os.writeKeyword("k") << k_ << token::END_STATEMENT << nl;
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

