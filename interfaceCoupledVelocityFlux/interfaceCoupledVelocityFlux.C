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
#include "interfaceCoupledVelocityValue.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCoupledVelocityFlux::
interfaceCoupledVelocityFlux
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    interfaceToInterfaceCoupleManager(p),
    kName_("k"),
    neighbourRegionName_(),
    neighbourPatchName_(),
    phiName_("phi"),
    rhoName_("rhoA"),
    nonOrthCorr_(false),
    secondOrder_(false)
{}


Foam::interfaceCoupledVelocityFlux::
interfaceCoupledVelocityFlux
(
    const interfaceCoupledVelocityFlux& icvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(icvf, p, iF, mapper),
    interfaceToInterfaceCoupleManager(p),
    kName_(icvf.kName_),
    neighbourRegionName_(icvf.neighbourRegionName_),
    neighbourPatchName_(icvf.neighbourPatchName_),
    phiName_(icvf.phiName_),
    rhoName_(icvf.rhoName_),
    nonOrthCorr_(icvf.nonOrthCorr_),
    secondOrder_(icvf.secondOrder_)
{}


Foam::interfaceCoupledVelocityFlux::
interfaceCoupledVelocityFlux
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    interfaceToInterfaceCoupleManager(p, dict),
    kName_(dict.lookupOrDefault<word>("k", word::null)),
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
    const interfaceCoupledVelocityFlux& icvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(icvf, iF),
    interfaceToInterfaceCoupleManager(icvf),
    kName_(icvf.kName_),
    neighbourRegionName_(icvf.neighbourRegionName_),
    neighbourPatchName_(icvf.neighbourPatchName_),
    phiName_(icvf.phiName_),
    rhoName_(icvf.rhoName_),
    nonOrthCorr_(icvf.nonOrthCorr_),
    secondOrder_(icvf.secondOrder_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceCoupledVelocityFlux::~interfaceCoupledVelocityFlux()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> interfaceCoupledVelocityFlux::velJump() const
{
    const areaVectorField& nf = rgInterface().aMesh().faceAreaNormals();                        

    // surface tension
    areaScalarField sigma = rgInterface().sigma();

    areaVectorField gradSsigma = fac::grad(sigma);

    // surface velocity terms
    const areaVectorField& Us = rgInterface().Us();
  
    areaScalarField divSU = fac::div(Us); 
    divSU.correctBoundaryConditions();

    areaTensorField gradSU = fac::grad(Us);

    // Remove component of gradient normal to surface (area)
    gradSU -= nf*(nf & gradSU);
    gradSU.correctBoundaryConditions();

    dimensionedScalar muFluidNbr
    (
        db().time().lookupObject<IOdictionary>("transportProperties")
        .subDict(nbrMesh().name()).lookup("mu")
    );
    
    dimensionedScalar muFluid
    (
        db().time().lookupObject<IOdictionary>("transportProperties")
        .subDict(patch().boundaryMesh().mesh().name()).lookup("mu")
    );
    
    return
    (
//        gradSsigma.internalField()
        (muFluidNbr.value() - muFluid.value())
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

    updateRegionInterface();

    // Calculate interpolated patch field
    vectorField fluxNbrToOwn(patch().size(), pTraits<vector>::zero);

    // Lookup neighbouring patch field
    const volVectorField& nbrField = nbrMesh().lookupObject<volVectorField>
        (
            //presume same field name as on this side
            this->dimensionedInternalField().name()
        );

    const fvPatchVectorField& nbrPatchField =
     (
         nbrPatch().patchField<volVectorField, vector>(nbrField)
     );
     
    tmp<vectorField> tnbrFlux = 
        refCast<const interfaceCoupledVelocityValue>
        (nbrPatchField).flux();
     
    const vectorField& nbrFlux = tnbrFlux();

    fluxNbrToOwn = interpolateFromNbrField<vector>(nbrFlux);

    fluxNbrToOwn *= -1.0;
    fluxNbrToOwn += velJump();

	// Enforce flux matching
    dimensionedScalar k
    (
        db().time().lookupObject<IOdictionary>("transportProperties")
        .subDict(patch().boundaryMesh().mesh().name()).lookup(kName_)
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
    // Calculate interpolated patch flux
    vectorField fluxNbrToOwn(patch().size(), pTraits<vector>::zero);

    // Lookup neighbouring patch field
    const volVectorField& nbrField = nbrMesh().lookupObject<volVectorField>
        (
            //presume same field name as on this side
            this->dimensionedInternalField().name()
        );
        
    const fvPatchVectorField& nbrPatchField =
        (
            nbrPatch().patchField<volVectorField, vector>(nbrField)
        );
     
    vectorField nbrFlux = 
        refCast<const interfaceCoupledVelocityValue>
        (nbrPatchField).flux();

    fluxNbrToOwn = interpolateFromNbrField<vector>(nbrFlux);

    // Calculate the maximum normalized residual
    const fvPatchVectorField& fown = *this;

    dimensionedScalar k
    (
        db().time().lookupObject<IOdictionary>("transportProperties")
        .subDict(patch().boundaryMesh().mesh().name()).lookup(kName_)
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
    os.writeKeyword("k") << kName_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourRegionName") << neighbourRegionName_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourPatchName") << neighbourPatchName_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("phiName") << phiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("rhoName") << rhoName_ << token::END_STATEMENT << nl;
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
    interfaceCoupledVelocityFlux
);

} // End namespace Foam

// ************************************************************************* //
