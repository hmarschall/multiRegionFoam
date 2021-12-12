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

#include "interfaceCoupledPressureValue.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "areaFields.H"
#include "mathematicalConstants.H"
#include "fvCFD.H"
#include "ggiInterpolation.H"

#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"
#include "wallFvPatch.H"
#include "fixedGradientFaPatchFields.H"

#include "wedgeFaPatch.H"
#include "wedgeFaPatchFields.H"
#include "slipFaPatchFields.H"

#include "interfaceCoupledVelocityValue.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::meshInterface& 
Foam::interfaceCoupledPressureValue::ale() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const objectRegistry& obr = mesh.objectRegistry::parent();

    if (!obr.foundObject<meshInterface>("aleProperties"))
    {
        FatalErrorIn("interfaceCoupledPressureValue::")
            << "meshInterface object not found but required."
            << abort(FatalError);
    }

    return obr.lookupObject<meshInterface>("aleProperties");
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCoupledPressureValue::
interfaceCoupledPressureValue
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    k_("k"),
    muFluidA_
    (
        patch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("muFluidA")
    ),
    muFluidB_
    (
        patch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("muFluidB")
    ),
    rhoFluidA_
    (
        patch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("rhoFluidA")
    ),
    rhoFluidB_
    (
        patch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("rhoFluidB")
    ),
    sigma_
    (
        patch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("cleanSurfaceTension")
    ),
    sigmaPtr_(NULL),
    g_
    (
        patch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("g")
    ),
    relax_(1.0),
    nonOrthCorr_(false),
    secondOrder_(false)
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
            this->patch().boundaryMesh().mesh()
            .lookupObject<IOdictionary>("surfaceProperties")
        )
    );
}


Foam::interfaceCoupledPressureValue::
interfaceCoupledPressureValue
(
    const interfaceCoupledPressureValue& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    k_(ptf.k_),
    muFluidA_(ptf.muFluidA_),
    muFluidB_(ptf.muFluidB_),
    rhoFluidA_(ptf.rhoFluidA_),
    rhoFluidB_(ptf.rhoFluidB_),
    sigma_(ptf.sigma_),
    sigmaPtr_(ptf.sigmaPtr_),
    g_(ptf.g_),
    relax_(ptf.relax_),
    nonOrthCorr_(ptf.nonOrthCorr_),
    secondOrder_(ptf.secondOrder_)
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
            this->patch().boundaryMesh().mesh()
            .lookupObject<IOdictionary>("surfaceProperties")
        )
    );
}


Foam::interfaceCoupledPressureValue::
interfaceCoupledPressureValue
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    k_(dict.lookupOrDefault<word>("k", word::null)),
    muFluidA_
    (
        patch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("muFluidA")
    ),
    muFluidB_
    (
        patch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("muFluidB")
    ),
    rhoFluidA_
    (
        patch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("rhoFluidA")
    ),
    rhoFluidB_
    (
        patch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("rhoFluidB")
    ),
    sigma_
    (
        patch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("cleanSurfaceTension")
    ),
    sigmaPtr_(NULL),
    g_
    (
        patch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("g")
    ),
    relax_(dict.lookupOrDefault<scalar>("relax",1.0)),
    nonOrthCorr_(false),
    secondOrder_(false)
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
            this->patch().boundaryMesh().mesh()
            .lookupObject<IOdictionary>("surfaceProperties")
        )
    );

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


Foam::interfaceCoupledPressureValue::
interfaceCoupledPressureValue
(
    const interfaceCoupledPressureValue& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wtcsf, iF),
    k_(wtcsf.k_),
    muFluidA_(wtcsf.muFluidA_),
    muFluidB_(wtcsf.muFluidB_),
    rhoFluidA_(wtcsf.rhoFluidA_),
    rhoFluidB_(wtcsf.rhoFluidB_),
    sigma_(wtcsf.sigma_),
    sigmaPtr_(wtcsf.sigmaPtr_),
    g_(wtcsf.g_),
    relax_(wtcsf.relax_),
    nonOrthCorr_(wtcsf.nonOrthCorr_),
    secondOrder_(wtcsf.secondOrder_)
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
            this->patch().boundaryMesh().mesh()
            .lookupObject<IOdictionary>("surfaceProperties")
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceCoupledPressureValue::~interfaceCoupledPressureValue()
{
    delete sigmaPtr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> interfaceCoupledPressureValue::valueJump() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const faMesh& aMesh = faMesh(mesh);

    const areaVectorField& nf = aMesh.faceAreaNormals();

    // Get interfacial curvature
    areaScalarField& K = interfacePtr_->K();
    interfacePtr_->updateK();


//    const fvMesh& meshA = ale().meshA();

//    if (!obr.foundObject<interfacialTransport>("surfaceProperties"))
//    {
//        FatalErrorIn("interfaceCoupledPressureValue::")
//            << "interfacialTransport object not found but required."
//            << abort(FatalError);
//    }

//    const interfacialTransport& interface =
//        meshA.lookupObject<interfacialTransport>("surfaceProperties");

    // surface velocity terms
    areaVectorField& Us = interfacePtr_->Us();
    interfacePtr_->updateUs();

    areaScalarField divUs
    (
        fac::div(Us)
    );
//    divUs.correctBoundaryConditions();

    // surface tension
    if 
    (
        mesh.foundObject<areaScalarField>("surfaceTension")
    )
    {
        // contaminated interface
        sigmaPtr_ = const_cast<areaScalarField*> 
            (&mesh.lookupObject<areaScalarField>("surfaceTension"));
    }
    else
    {
        // clean interface
        sigmaPtr_ = new areaScalarField
        (
            IOobject
            (
                "sigma",
                this->db().time().timeName(),
                aMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh,
            sigma_,
            zeroGradientFaPatchVectorField::typeName
        );
    }

    areaScalarField& sigma = *sigmaPtr_;

    // gravity term
    vector pRefPoint(mesh.solutionDict().subDict("PISO").lookup("pRefPoint"));

    return
    (
        2.0*(muFluidB_.value() - muFluidA_.value())*divUs.internalField()
      - sigma.internalField()*K.internalField()
      + (rhoFluidB_.value() - rhoFluidA_.value())
        *(
            (
                mesh.C().boundaryField()[this->patch().index()] 
              - pRefPoint
            ) & g_.value()
        )
   );
}

//- Update the patch field coefficients
void Foam::interfaceCoupledPressureValue::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Calculate interpolated patch field
    scalarField fieldNbrToOwn(patch().size(), 0);

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


    //- Lookup neighbouring patch field
    const volScalarField& nbrField = ale().meshB().lookupObject<volScalarField>
        (
            //presume same field name as on this side
            this->dimensionedInternalField().name()
        );

    const fvPatchScalarField& nbrPatchField =
    (
        ale().patchB().patchField<volScalarField, scalar>(nbrField)
    );

    //- Interpolate pressure face values B-to-A
    ale().transferFacesBToA(nbrPatchField, fieldNbrToOwn);

    // Add interfacial pressure jump
    fieldNbrToOwn += valueJump();

	operator==
    (
        *this
      + relax_*
        (
            fieldNbrToOwn
          - *this
//          - sqr(interfacePtr_->massFluxCorr())
//            *(
//                 1.0/interfacePtr_->rhoB().value()
//               - 1.0/interfacePtr_->rhoA().value()
//             )
        )
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


//- Return patch-normal gradient
Foam::tmp<Foam::Field<Foam::scalar> >
Foam::interfaceCoupledPressureValue::snGrad() const
{
    word pName = this->dimensionedInternalField().name();

    const fvPatchField<vector>& gradp =
        patch().lookupPatchField<volVectorField, vector>
        (
            "grad(" + pName + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = delta - n*(n&delta);

    scalarField dpP(this->patch().size(), 0);
    if (nonOrthCorr_)
    {
        dpP = (k&gradp.patchInternalField());
    }

    if (secondOrder_)
    {
        scalarField nGradpP = (n&gradp.patchInternalField());

        return
            2
           *(
                *this
              - (this->patchInternalField() + dpP)
            )*this->patch().deltaCoeffs()
          - nGradpP;
    }

    return
    (
        *this
      - (patchInternalField() + (k&gradp.patchInternalField()))
    )*this->patch().deltaCoeffs();
}

//- Return the patch flux
Foam::tmp<Foam::scalarField>
Foam::interfaceCoupledPressureValue::flux() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    dimensionedScalar k
    (
      mesh.lookupObject<IOdictionary>
      ("surfaceProperties").lookup(k_)
    );

    return (this->snGrad()/k.value());
}

Foam::tmp<Foam::Field<Foam::scalar> > 
Foam::interfaceCoupledPressureValue::gradientBoundaryCoeffs() const
{
    word pName = this->dimensionedInternalField().name();

    const fvPatchField<vector>& gradp =
        patch().lookupPatchField<volVectorField, vector>
        (
            "grad(" + pName + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = delta - n*(n&delta);

    scalarField dpP(this->patch().size(), 0);
    if (nonOrthCorr_)
    {
        dpP = (k&gradp.patchInternalField());
    }

    if (secondOrder_)
    {
        scalarField nGradpP = (n&gradp.patchInternalField());

        return
            this->patch().deltaCoeffs()
           *(
                2*(*this - dpP) 
              - this->patchInternalField()
            )
          - nGradpP;
    }

    return this->patch().deltaCoeffs()*(*this - dpP);
}

//- Return the maximum normalized coupled patch residual
Foam::scalarField Foam::interfaceCoupledPressureValue::residual() const
{
    // Calculate interpolated patch field
    const polyPatch& ownPatch = patch().patch();
    const polyPatch& nbrPatch = coupleManagerPtr_->neighbourPatch().patch();

    const GGIInterpolation<polyPatch, polyPatch>& interpolator = 
        coupleManagerPtr_->ggiInterpolator(nbrPatch,ownPatch);

    const fvPatchScalarField& nbrPatchField =
         coupleManagerPtr_->neighbourPatchField<scalar>();

    scalarField fieldNbrToOwn = interpolator.masterToSlave
    (
        nbrPatchField
    );

	// Calculate the maximal normalized residual
	const scalarField& fown = *this;

    scalarField K
    (
        valueJump()/(fieldNbrToOwn + SMALL)
//        valueJump()/(*this + SMALL)
    );

    const scalarField& residual =
        (
            (mag(fown - fieldNbrToOwn) - mag(valueJump()))/
            min(mag(fown),mag(K*fieldNbrToOwn))
        );

//    const scalar& residual =
//        gMax
//        (
//            mag(fown - K*fieldNbrToOwn)/
//            max(min(gMax(fown),gMax(K*fieldNbrToOwn)), SMALL)
//        );

	return residual;
}


//- Write
void Foam::interfaceCoupledPressureValue::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    coupleManagerPtr_->writeEntries(os);
    os.writeKeyword("k") << k_ << token::END_STATEMENT << nl;
    os.writeKeyword("relax") << relax_ << token::END_STATEMENT << nl;
    os.writeKeyword("nonOrthCorr") << nonOrthCorr_ << token::END_STATEMENT << nl;
    os.writeKeyword("secondOrder") << secondOrder_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    interfaceCoupledPressureValue
);

} // End namespace Foam

// ************************************************************************* //
