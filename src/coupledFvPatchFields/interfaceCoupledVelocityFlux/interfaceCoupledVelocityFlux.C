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

#include "interfaceCoupledVelocityFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"
#include "interfaceCoupledVelocityValue.H"
#include "ggiInterpolation.H"
#include "patchToPatchInterpolation.H"

#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"
#include "wallFvPatch.H"

#include "wedgeFaPatch.H"
#include "wedgeFaPatchFields.H"
#include "slipFaPatchFields.H"

#include "fixedGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::meshInterface& 
Foam::interfaceCoupledVelocityFlux::ale() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const objectRegistry& obr = mesh.objectRegistry::parent();

    if (!obr.foundObject<meshInterface>("aleProperties"))
    {
        FatalErrorIn("interfaceCoupledVelocityFlux::")
            << "meshInterface object not found but required."
            << abort(FatalError);
    }

    return obr.lookupObject<meshInterface>("aleProperties");
}

const Foam::interfacialTransport& 
Foam::interfaceCoupledVelocityFlux::interface() const
{
    const fvMesh& mesh = ale().meshA();

    if (!mesh.foundObject<interfacialTransport>("interfaceTransportProperties"))
    {
        FatalErrorIn("interfaceCoupledVelocityFlux::")
            << "interfacialTransport object not found but required."
            << abort(FatalError);
    }

    return
    (
        mesh.lookupObject<interfacialTransport>
        ("interfaceTransportProperties")
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCoupledVelocityFlux::
interfaceCoupledVelocityFlux
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
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
    sigma_
    (
        patch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("cleanSurfaceTension")
    ),
    sigmaPtr_(NULL),
    phiName_("phi"),
    rhoName_("rhoA"),
    nonOrthCorr_(false),
    secondOrder_(false),
    curTimeIndex_( dimensionedInternalField().mesh().time().timeIndex() ),
    Fc_( p.patch().size(), vector::zero ),
    oldFc_( p.patch().size(), vector::zero ),
    oldoldFc_( p.patch().size(), vector::zero )
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


Foam::interfaceCoupledVelocityFlux::
interfaceCoupledVelocityFlux
(
    const interfaceCoupledVelocityFlux& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(ptf, p, iF, mapper),
    k_(ptf.k_),
    muFluidA_(ptf.muFluidA_),
    muFluidB_(ptf.muFluidB_),
    sigma_(ptf.sigma_),
    sigmaPtr_(NULL),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    nonOrthCorr_(ptf.nonOrthCorr_),
    secondOrder_(ptf.secondOrder_),
    curTimeIndex_(ptf.curTimeIndex_),
    Fc_(p.patch().size(), vector::zero),
    oldFc_(p.patch().size(), vector::zero),
    oldoldFc_(p.patch().size(), vector::zero)
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


Foam::interfaceCoupledVelocityFlux::
interfaceCoupledVelocityFlux
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
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
    sigma_
    (
        patch().boundaryMesh().mesh()
        .lookupObject<IOdictionary>("surfaceProperties")
        .lookup("cleanSurfaceTension")
    ),
    sigmaPtr_(NULL),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    nonOrthCorr_(false),
    secondOrder_(false),
    curTimeIndex_(dimensionedInternalField().mesh().time().timeIndex()),
    Fc_(p.patch().size(), vector::zero),
    oldFc_(p.patch().size(), vector::zero),
    oldoldFc_(p.patch().size(), vector::zero)
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

    fvPatchVectorField::operator=(patchInternalField());

    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
        fvPatchField<vector>::operator=(this->patchInternalField()());
    }
    else
    {
        evaluate();
    }
}


Foam::interfaceCoupledVelocityFlux::
interfaceCoupledVelocityFlux
(
    const interfaceCoupledVelocityFlux& whftcsf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(whftcsf, iF),
//    coupleManager_(whftcsf.coupleManager_),
    k_(whftcsf.k_),
    muFluidA_(whftcsf.muFluidA_),
    muFluidB_(whftcsf.muFluidB_),
    sigma_(whftcsf.sigma_),
    sigmaPtr_(NULL),
    phiName_(whftcsf.phiName_),
    rhoName_(whftcsf.rhoName_),
    nonOrthCorr_(whftcsf.nonOrthCorr_),
    secondOrder_(whftcsf.secondOrder_),
    curTimeIndex_(whftcsf.curTimeIndex_),
    Fc_(whftcsf.oldFc_),
    oldFc_(whftcsf.oldFc_),
    oldoldFc_(whftcsf.oldoldFc_)
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
            this->patch().name(),
            this->patch().boundaryMesh().mesh()
            .lookupObject<IOdictionary>("surfaceProperties")
        )
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceCoupledVelocityFlux::~interfaceCoupledVelocityFlux()
{
    delete sigmaPtr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> interfaceCoupledVelocityFlux::velJump() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const faMesh& aMesh = interfacePtr_->aMesh();

    const areaVectorField& nf = aMesh.faceAreaNormals();

    // Get interfacial curvature
//    areaScalarField& K = interfacePtr_->K();
//    interfacePtr_->updateK();

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

    areaVectorField gradSsigma =
    (
        fac::grad(sigma)
    );

    // surface velocity terms
    areaVectorField& Us = interfacePtr_->Us();
    interfacePtr_->updateUs();

    areaScalarField divSU =
    (
        fac::div(Us)
    );

    areaTensorField gradSU =
    (
        fac::grad(Us)
    );

    // Remove component of gradient normal to surface (area)
    gradSU -= nf*(nf & gradSU);
    gradSU.correctBoundaryConditions();

    return
    (
//        gradSsigma.internalField()
        (muFluidB_.value() - muFluidA_.value())
        *(
            divSU*nf
          + (gradSU&nf)
        )().internalField()
    );
}


//- Update the patch field coefficients
void Foam::interfaceCoupledVelocityFlux::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // Calculate interpolated patch field
    vectorField fluxNbrToOwn(patch().size(), pTraits<vector>::zero);
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
    const volVectorField& nbrField = ale().meshB().lookupObject<volVectorField>
        (
            //presume same field name as on this side
            this->dimensionedInternalField().name()
        );

    const fvPatchVectorField& nbrPatchField =
    (
        ale().patchB().patchField<volVectorField, vector>(nbrField)
    );

    //- Interpolate flux face values B-to-A
    tmp<vectorField> tnbrFlux = 
        refCast<const interfaceCoupledVelocityValue>
        (nbrPatchField).flux();
    const vectorField& nbrFlux = tnbrFlux();

    ale().transferFacesBToA
    (
        nbrFlux,
        fluxNbrToOwn
    );

    fluxNbrToOwn *= -1.0;
    fluxNbrToOwn += velJump();

	// Enforce flux matching
    dimensionedScalar k
    (
      mesh.lookupObject<IOdictionary>
      ("surfaceProperties").lookup(k_)
    );

    gradient() = fluxNbrToOwn/k.value();

    fixedGradientFvPatchVectorField::updateCoeffs();
}

void Foam::interfaceCoupledVelocityFlux::evaluate
(
    const Pstream::commsTypes
)
{
    if (!updated())
    {
        updateCoeffs();
    }

    fixedGradientFvPatchVectorField::evaluate();

    // Evaluate tangential component of velocity
    // using (second order discretisation) and nonorthogonal correction
    vectorField dUP(this->patch().size(), vector::zero);

    {
        vectorField n = patch().nf();
        vectorField delta = patch().delta();
        vectorField k = delta - n*(n&delta);

        word UName = this->dimensionedInternalField().name();

        const fvPatchField<tensor>& gradU =
            patch().lookupPatchField<volTensorField, tensor>
            (
                "grad(" + UName + ")"
            );
        if (nonOrthCorr_)
        {
            dUP = (k&gradU.patchInternalField());
        }

//        surfaceScalarField& phiField = const_cast<surfaceScalarField&>
//        (
//            this->db().lookupObject<surfaceScalarField>(phiName_)
//        );

        fvsPatchField<scalar>& phiFieldp = const_cast<fvsPatchField<scalar>& >
        (
            this->db().lookupObject<surfaceScalarField>(phiName_)
            .boundaryField()[this->patch().index()]
        );

        if (secondOrder_)
        {
            vectorField nGradUP = (n&gradU.patchInternalField());

            fvPatchField<vector>::operator=
            (
                this->patchInternalField() + dUP
              + 0.5*(gradient() + nGradUP)
                /this->patch().deltaCoeffs()
//              - interfacePtr_->massFluxCorr()
//                *(
//                    1./interfacePtr_->rhoB()
//                  - 1./interfacePtr_->rhoA()
//                ).value()*patch().nf()
            );

            phiFieldp = 
            (
                (
                    (
                        this->patchInternalField() + dUP
                      + 0.5*(gradient() + nGradUP)
                        /this->patch().deltaCoeffs()
                    ) 
                    & patch().Sf()
                )
            );
        }
        else
        {
            fvPatchField<vector>::operator=
            (
                this->patchInternalField() + dUP
              + gradient()/this->patch().deltaCoeffs()
//              - interfacePtr_->massFluxCorr()
//                *(
//                    1./interfacePtr_->rhoB()
//                  - 1./interfacePtr_->rhoA()
//                ).value()*patch().nf()
            );

            phiFieldp = 
            (
                (
                    (
                        this->patchInternalField() + dUP
                      + gradient()/this->patch().deltaCoeffs()
                    ) 
                    & patch().Sf()
                )
            );
        }
    }

    fvPatchField<vector>::evaluate();
}

//- Return the maximum normalized coupled patch residual
Foam::scalarField Foam::interfaceCoupledVelocityFlux::residual() const
{

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // Calculate interpolated patch flux
    const polyPatch& ownPatch = patch().patch();
    const polyPatch& nbrPatch = coupleManagerPtr_->neighbourPatch().patch();

    const fvPatchField<vector>& nbrPatchField =
        coupleManagerPtr_->neighbourPatchField<vector>();

    const GGIInterpolation<polyPatch, polyPatch>& interpolator =
        coupleManagerPtr_->ggiInterpolator(nbrPatch,ownPatch);

    vectorField fluxNbrToOwn =
        interpolator.masterToSlave
        (
            refCast<const interfaceCoupledVelocityValue>
            (nbrPatchField).flux()
        );

    // Calculate the maximum normalized residual
    const fvPatchVectorField& fown = *this;

    dimensionedScalar k
    (
        mesh.lookupObject<IOdictionary>
        ("surfaceProperties").lookup(k_)
    );

    vectorField fluxOwn = k.value()*fown.snGrad();

    const scalarField& residual = 
        (mag(mag(fluxOwn) - mag(fluxNbrToOwn)) - mag(velJump()))/
        min(mag(fluxOwn), mag(fluxNbrToOwn));
//        (max(min(gMax(mag(fluxOwn)), gMax(mag(fluxNbrToOwn))), SMALL));

    //scalarField residual =
  	//    mag(mag(fluxOwn) - mag(fluxNbrToOwn))/
    //          (max(min(gMax(mag(fluxOwn)),gMax(mag(fluxNbrToOwn))), SMALL ));

    return residual;
}

//- Write
void Foam::interfaceCoupledVelocityFlux::write
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
    interfaceCoupledVelocityFlux
);

} // End namespace Foam

// ************************************************************************* //
