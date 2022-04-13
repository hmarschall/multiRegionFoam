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
#include "interfaceCoupledPressureValue.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCoupledPressureFlux::
interfaceCoupledPressureFlux
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    interfaceToInterfaceCoupleManager(p),
    kName_("k"),
    neighbourRegionName_(),
    neighbourPatchName_()
{}


Foam::interfaceCoupledPressureFlux::
interfaceCoupledPressureFlux
(
    const interfaceCoupledPressureFlux& icpf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(icpf, p, iF, mapper),
    interfaceToInterfaceCoupleManager(p),
    kName_(icpf.kName_),
    neighbourRegionName_(icpf.neighbourRegionName_),
    neighbourPatchName_(icpf.neighbourPatchName_)
{}


Foam::interfaceCoupledPressureFlux::
interfaceCoupledPressureFlux
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    interfaceToInterfaceCoupleManager(p, dict),
    kName_(dict.lookupOrDefault<word>("k", word::null)),
    neighbourRegionName_
    (
        dict.lookupOrDefault<word>("neighbourRegionName", word::null)
    ),
    neighbourPatchName_
    (
        dict.lookupOrDefault<word>("neighbourPatchName", word::null)
    )
{
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
    const interfaceCoupledPressureFlux& icpf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(icpf, iF),
    interfaceToInterfaceCoupleManager(icpf),
    kName_(icpf.kName_),
    neighbourRegionName_(icpf.neighbourRegionName_),
    neighbourPatchName_(icpf.neighbourPatchName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> interfaceCoupledPressureFlux::nGradJump() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const volVectorField& U = 
        mesh.objectRegistry::lookupObject<volVectorField>("U");

    const volVectorField& nbrU = 
        nbrMesh().lookupObject<volVectorField>("U");

    const surfaceScalarField& phi = 
        mesh.objectRegistry::lookupObject<surfaceScalarField>("phi");

    const surfaceScalarField& nbrPhi = 
        nbrMesh().lookupObject<surfaceScalarField>("phi");

    vectorField nbrLaplacianU = interpolateFromNbrField<vector> 
        (
            fvc::laplacian(nbrU)() 
            .boundaryField()[nbrPatch().index()] 
            .patchInternalField()
        );

    vectorField nbrDDtU = interpolateFromNbrField<vector> 
        (
            fvc::DDt(nbrPhi, nbrU)()
            .boundaryField()[nbrPatch().index()] 
            .patchInternalField()
        );
    
    dimensionedScalar muFluidNbr
    (
        db().time().lookupObject<IOdictionary>("transportProperties")
        .subDict(nbrMesh().name()).lookup("mu")
    );
    
    dimensionedScalar muFluid 
    (
        db().time().lookupObject<IOdictionary>("transportProperties")
        .subDict(mesh.name()).lookup("mu")
    );
    
    dimensionedScalar rhoFluidNbr 
    (
        db().time().lookupObject<IOdictionary>("transportProperties")
        .subDict(nbrMesh().name()).lookup("rho")
    );
    
    dimensionedScalar rhoFluid 
    (
        db().time().lookupObject<IOdictionary>("transportProperties")
        .subDict(mesh.name()).lookup("rho")
    );
    
    return 
    (
        (
            muFluid.value()/rhoFluid.value()
            *(
                fvc::laplacian(U)()
                .boundaryField()[this->patch().index()]
                .patchInternalField()
            )
          - muFluidNbr.value()/rhoFluidNbr.value()*
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

    // Lookup neighbouring patch field
    const volScalarField& nbrField = nbrMesh().lookupObject<volScalarField>
        (
            //presume same field name as on this side
            this->dimensionedInternalField().name()
        );

    // Interpolate flux face values from neighbour patch
    tmp<scalarField> tnbrFlux = 
        refCast<const interfaceCoupledPressureValue>
        (nbrPatch().patchField<volScalarField, scalar>(nbrField)).flux(); 

    const scalarField& nbrFlux = tnbrFlux();

    fluxNbrToOwn = interpolateFromNbrField<scalar>(nbrFlux);

	// Enforce flux matching
    fluxNbrToOwn *= -1.0; 
    fluxNbrToOwn += nGradJump();

    dimensionedScalar k
    (
        db().time().lookupObject<IOdictionary>("transportProperties")
        .subDict(patch().boundaryMesh().mesh().name()).lookup(kName_)
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

    word fieldName = dimensionedInternalField().name();

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
    // Calculate interpolated patch flux
    scalarField fluxNbrToOwn(patch().size(), 0);
    
    // Lookup neighbouring patch field
    const volScalarField& nbrField = nbrMesh().lookupObject<volScalarField>
    (
        //presume same field name as on this side
        this->dimensionedInternalField().name()
    );

    tmp<scalarField> tnbrFlux = 
        refCast<const interfaceCoupledPressureValue>
        (nbrPatch().patchField<volScalarField, scalar>(nbrField)).flux();  
    const scalarField& nbrFlux = tnbrFlux();

    fluxNbrToOwn = interpolateFromNbrField<scalar>(nbrFlux);

    // Calculate the maximum normalized residual
    const fvPatchScalarField& fown = *this;
    
    dimensionedScalar k
    (
        db().time().lookupObject<IOdictionary>("transportProperties")
        .subDict(patch().boundaryMesh().mesh().name()).lookup(kName_)
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
    os.writeKeyword("k") << kName_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourRegionName") << neighbourRegionName_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourPatchName") << neighbourPatchName_ 
        << token::END_STATEMENT << nl;    
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
