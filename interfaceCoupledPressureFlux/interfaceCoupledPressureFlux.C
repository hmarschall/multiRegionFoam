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

#include "interfaceCoupledPressureFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"
#include "interfaceCoupledPressureValue.H"
#include "ggiInterpolation.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::meshInterface& 
Foam::interfaceCoupledPressureFlux::ale() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const objectRegistry& obr = mesh.objectRegistry::parent();

    if (!obr.foundObject<meshInterface>("aleProperties"))
    {
        FatalErrorIn("interfaceCoupledPressureFlux::")
            << "meshInterface object not found but required."
            << abort(FatalError);
    }

    return obr.lookupObject<meshInterface>("aleProperties");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCoupledPressureFlux::
interfaceCoupledPressureFlux
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    k_("k"),
    muFluidA_
    (
        dimensionedScalar("muFluidA", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0)
    ),
    muFluidB_
    (
        dimensionedScalar("muFluidB", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0)
    ),
    rhoFluidA_
    (
        dimensionedScalar("rhoFluidA", dimMass/dimVolume, 0)
    ),
    rhoFluidB_
    (
        dimensionedScalar("rhoFluidB", dimMass/dimVolume, 0)
    ),
    g_
    (
        dimensionedVector("g", dimensionSet(0, 1, -2, 0, 0, 0, 0), vector::zero)
    )
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
//            coupleManagerPtr_->neighbourPatch().name(),
            this->patch().name(),
            coupleManagerPtr_->neighbourPatch().boundaryMesh().mesh()
            .lookupObject<IOdictionary>("surfaceProperties")
        )
    );

    muFluidA_ = dimensionedScalar
    (
        coupleManagerPtr_->neighbourPatch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("muFluidA")
    );

    muFluidB_ = dimensionedScalar
    (
        coupleManagerPtr_->neighbourPatch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("muFluidB")
    );

    rhoFluidA_ = dimensionedScalar
    (
        coupleManagerPtr_->neighbourPatch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("rhoFluidA")
    );

    rhoFluidB_ = dimensionedScalar
    (
        coupleManagerPtr_->neighbourPatch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("rhoFluidB")
    );

    g_ = dimensionedVector
    (
        coupleManagerPtr_->neighbourPatch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("g")
    );
}


Foam::interfaceCoupledPressureFlux::
interfaceCoupledPressureFlux
(
    const interfaceCoupledPressureFlux& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    k_(ptf.k_),
    muFluidA_(ptf.muFluidA_),
    muFluidB_(ptf.muFluidB_),
    rhoFluidA_(ptf.rhoFluidA_),
    rhoFluidB_(ptf.rhoFluidB_),
    g_(ptf.g_)
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
//            coupleManagerPtr_->neighbourPatch().name(),
            this->patch().name(),
            coupleManagerPtr_->neighbourPatch().boundaryMesh().mesh()
            .lookupObject<IOdictionary>("surfaceProperties")
        )
    );
}


Foam::interfaceCoupledPressureFlux::
interfaceCoupledPressureFlux
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    k_(dict.lookupOrDefault<word>("k", word::null)),
    muFluidA_
    (
        dimensionedScalar("muFluidA", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0)
    ),
    muFluidB_
    (
        dimensionedScalar("muFluidB", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0)
    ),
    rhoFluidA_
    (
        dimensionedScalar("rhoFluidA", dimMass/dimVolume, 0)
    ),
    rhoFluidB_
    (
        dimensionedScalar("rhoFluidB", dimMass/dimVolume, 0)
    ),
    g_
    (
        dimensionedVector("g", dimensionSet(0, 1, -2, 0, 0, 0, 0), vector::zero)
    )
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
//            coupleManagerPtr_->neighbourPatch().name(),
            this->patch().name(),
            coupleManagerPtr_->neighbourPatch().boundaryMesh().mesh()
            .lookupObject<IOdictionary>("surfaceProperties")
        )
    );

    muFluidA_ = dimensionedScalar
    (
        coupleManagerPtr_->neighbourPatch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("muFluidA")
    );

    muFluidB_ = dimensionedScalar
    (
        coupleManagerPtr_->neighbourPatch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("muFluidB")
    );

    rhoFluidA_ = dimensionedScalar
    (
        coupleManagerPtr_->neighbourPatch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("rhoFluidA")
    );

    rhoFluidB_ = dimensionedScalar
    (
        coupleManagerPtr_->neighbourPatch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("rhoFluidB")
    );

    g_ = dimensionedVector
    (
        coupleManagerPtr_->neighbourPatch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("g")
    );

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
        fvPatchField<scalar>::operator=(this->patchInternalField()());
    }
    else
    {
        evaluate();
    }
}


Foam::interfaceCoupledPressureFlux::
interfaceCoupledPressureFlux
(
    const interfaceCoupledPressureFlux& whftcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(whftcsf, iF),
    k_(whftcsf.k_),
    muFluidA_(whftcsf.muFluidA_),
    muFluidB_(whftcsf.muFluidB_),
    rhoFluidA_(whftcsf.rhoFluidA_),
    rhoFluidB_(whftcsf.rhoFluidB_),
    g_(whftcsf.g_)
{
    coupleManagerPtr_.reset
    (
        new patchCoupleManager
        (
            whftcsf.coupleManagerPtr_
        )
    );

    interfacePtr_.reset
    (
        new interfacialTransport
        (
//            word::null,
            this->db().lookupObject<volVectorField>("U"),
            this->db().lookupObject<surfaceScalarField>("phi"),
//            coupleManagerPtr_->neighbourPatch().name(),
            this->patch().name(),
            coupleManagerPtr_->neighbourPatch().boundaryMesh().mesh()
            .lookupObject<IOdictionary>("surfaceProperties")
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> interfaceCoupledPressureFlux::nGradJump() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
//    const objectRegistry& obr = mesh.objectRegistry::parent();

//    if (!obr.foundObject<meshInterface>("aleProperties"))
//    {
//        FatalErrorIn("interfaceCoupledPressureValue::updateCoeffs(...)")
//            << "meshInterface object not found but required."
//            << abort(FatalError);
//    }

//    const meshInterface& ale =
//        obr.lookupObject<meshInterface>("aleProperties");

    // Get interface normal
    const volVectorField& U =
        mesh.objectRegistry::lookupObject<volVectorField>("U");

    const volVectorField& nbrU = 
        ale().meshA().lookupObject<volVectorField>("U");

    const surfaceScalarField& phi =
        mesh.objectRegistry::lookupObject<surfaceScalarField>("phi");

    const surfaceScalarField& nbrPhi = 
        ale().meshA().lookupObject<surfaceScalarField>("phi");

    vectorField laplacianUA = fvc::laplacian(nbrU)()
        .boundaryField()[ale().patchAID()]
        .patchInternalField();

    vectorField DDtUA = fvc::DDt(nbrPhi, nbrU)()
        .boundaryField()[ale().patchAID()]
        .patchInternalField();

    vectorField nbrLaplacianU =
        ale().transferFacesFromA
        (
            laplacianUA
        );

    vectorField nbrDDtU =
        ale().transferFacesFromA
        (
            DDtUA
        );

    return 
    (
        (
            muFluidB_.value()/rhoFluidB_.value()
            *(
                fvc::laplacian(U)()
                .boundaryField()[this->patch().index()]
                .patchInternalField()
            )
          - muFluidA_.value()/rhoFluidA_.value()*
            nbrLaplacianU
        )
      + (
            nbrDDtU
          - (
                fvc::DDt(phi,U)()
                .boundaryField()[this->patch().index()]
                .patchInternalField()
            )
        )
    ) & patch().nf();
}

//- Update the patch field coefficients
void Foam::interfaceCoupledPressureFlux::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Calculate interpolated patch field
    scalarField fluxNbrToOwn(patch().size(), 0);

    //- Lookup neighbouring patch field
    const volScalarField& nbrField = ale().meshA().lookupObject<volScalarField>
        (
            //presume same field name as on this side
            this->dimensionedInternalField().name()
        );

    const fvPatchScalarField& nbrPatchField =
    (
        ale().patch().patchField<volScalarField, scalar>(nbrField)
    );

    //- Interpolate flux face values A-to-B
    tmp<scalarField> tnbrFlux = 
        refCast<const interfaceCoupledPressureValue>
        (nbrPatchField).flux();
    const scalarField& nbrFlux = tnbrFlux();

    ale().transferFacesAToB
    (
        nbrFlux,
        fluxNbrToOwn
    );

	// Enforce flux matching
    fluxNbrToOwn *= -1.0;
    fluxNbrToOwn += nGradJump();

    dimensionedScalar k
    (
      ale().meshA().lookupObject<IOdictionary>
      ("surfaceProperties").lookup(k_)
    );

    gradient() = fluxNbrToOwn*k.value();

    fixedGradientFvPatchScalarField::updateCoeffs();
}

void interfaceCoupledPressureFlux::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();
    word fieldName = dimensionedInternalField().name();

    const volScalarField& p =
        mesh.lookupObject<volScalarField>(fieldName);

    const fvPatchField<vector>& gradP =
        patch().lookupPatchField<volVectorField, vector>("grad(p)");

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = delta - n*(n&delta);

    scalarField dPP = (k&gradP.patchInternalField());

    scalarField nGradPb(this->patch().size(), 0);

    scalarField gradPp = (n&gradP.patchInternalField());

    //TODO: introduce and use secondOrder_ and nonOrthCorr_

    // second order with non-orthogonality correction
    Field<scalar>::operator=
    (
        this->patchInternalField() + dPP
//      + gradient()/this->patch().deltaCoeffs()
      + 0.5*(gradient() + gradPp)/this->patch().deltaCoeffs()
//      + nGradPb/this->patch().deltaCoeffs()
//      - sqr(interfacePtr_->massFluxCorr())
//        *(
//             1.0/interfacePtr_->rhoB().value()
//           - 1.0/interfacePtr_->rhoA().value()
//         )
    );

    fvPatchField<scalar>::evaluate();
}

//- Return the maximum normalized coupled patch residual
Foam::scalarField Foam::interfaceCoupledPressureFlux::residual() const
{
    const fvMesh& nbrMesh = coupleManagerPtr_->neighbourPatch().boundaryMesh().mesh();

    // Calculate interpolated patch flux
    const polyPatch& ownPatch = patch().patch();
    const polyPatch& nbrPatch = coupleManagerPtr_->neighbourPatch().patch();

    const GGIInterpolation<polyPatch, polyPatch>& interpolator = 
        coupleManagerPtr_->ggiInterpolator(nbrPatch,ownPatch);

    const fvPatchScalarField& nbrPatchField =
         coupleManagerPtr_->neighbourPatchField<scalar>();

    scalarField fluxNbrToOwn =
        interpolator.masterToSlave
        (
            refCast<const interfaceCoupledPressureValue>
            (nbrPatchField).flux()
        );

    // Calculate the maximum normalized residual
    const fvPatchScalarField& fown = *this;

    dimensionedScalar k
    (
      nbrMesh.lookupObject<IOdictionary>
      ("surfaceProperties").lookup(k_)
    );

    scalarField fluxOwn = 1.0/k.value()*fown.snGrad();

    const scalarField& residualField =
	    mag(nGradJump())/
        (max(min(gMax(mag(fluxOwn)), gMax(mag(fluxNbrToOwn))), SMALL));

    return residualField;
}


//- Write
void Foam::interfaceCoupledPressureFlux::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    coupleManagerPtr_->writeEntries(os);
    os.writeKeyword("k") << k_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    interfaceCoupledPressureFlux
);

} // End namespace Foam

// ************************************************************************* //
