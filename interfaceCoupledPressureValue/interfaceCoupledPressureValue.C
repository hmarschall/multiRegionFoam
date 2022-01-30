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
#include "interfaceCoupledVelocityValue.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCoupledPressureValue::
interfaceCoupledPressureValue
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    interfaceToInterfaceCoupleManager(p),
    kName_("k"),
    neighbourRegionName_(),
    neighbourPatchName_(), 
    relax_(1.0),
    nonOrthCorr_(false),
    secondOrder_(false)
{}


Foam::interfaceCoupledPressureValue::
interfaceCoupledPressureValue
(
    const interfaceCoupledPressureValue& icpv,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(icpv, p, iF, mapper),
    interfaceToInterfaceCoupleManager(p),
    kName_(icpv.kName_),
    neighbourRegionName_(icpv.neighbourRegionName_),
    neighbourPatchName_(icpv.neighbourPatchName_),
    relax_(icpv.relax_),
    nonOrthCorr_(icpv.nonOrthCorr_),
    secondOrder_(icpv.secondOrder_)
{}


Foam::interfaceCoupledPressureValue::
interfaceCoupledPressureValue
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
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
    relax_(dict.lookupOrDefault<scalar>("relax",1.0)),
    nonOrthCorr_(false),
    secondOrder_(false)
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


Foam::interfaceCoupledPressureValue::
interfaceCoupledPressureValue
(
    const interfaceCoupledPressureValue& icpv,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(icpv, iF),
    interfaceToInterfaceCoupleManager(icpv),
    kName_(icpv.kName_),
    neighbourRegionName_(icpv.neighbourRegionName_),
    neighbourPatchName_(icpv.neighbourPatchName_),
    relax_(icpv.relax_),
    nonOrthCorr_(icpv.nonOrthCorr_),
    secondOrder_(icpv.secondOrder_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceCoupledPressureValue::~interfaceCoupledPressureValue()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> interfaceCoupledPressureValue::valueJump() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const faMesh& aMesh = faMesh(mesh);

    // interfacial curvature
    const areaScalarField& K = rgInterface().K();
   
    // surface velocity terms
    const areaVectorField& Us = rgInterface().Us();

    areaScalarField divUs(fac::div(Us));
    
//    divUs.correctBoundaryConditions();

    // surface tension 
    areaScalarField sigma = rgInterface().sigma();  

    // gravity term
    vector pRefPoint(mesh.solutionDict().subDict("PISO").lookup("pRefPoint"));
    
    dimensionedVector g (rgInterface().gravitationalProperties().lookup("g"));

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
        2.0*(muFluidNbr.value() - muFluid.value())*divUs.internalField()
      - sigma.internalField()*K.internalField()
      + (rhoFluidNbr.value() - rhoFluid.value())
        *(
            (
                mesh.C().boundaryField()[this->patch().index()] 
              - pRefPoint
            ) & g.value()
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

    updateRegionInterface();
    
    // Calculate interpolated patch field
    scalarField fieldNbrToOwn(patch().size(), 0);

    // Lookup neighbouring patch field
    const volScalarField& nbrField = nbrMesh().lookupObject<volScalarField>
        (
            //presume same field name as on this side
            this->dimensionedInternalField().name()
        );

    fieldNbrToOwn = interpolateFromNbrField<scalar>
        (
            nbrPatch().patchField<volScalarField, scalar>(nbrField)
        );

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
    dimensionedScalar k
    (
        db().time().lookupObject<IOdictionary>("transportProperties")
        .subDict(patch().boundaryMesh().mesh().name()).lookup(kName_)
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
    scalarField fieldNbrToOwn(patch().size(), 0);
    
    // Lookup neighbouring patch field
    const volScalarField& nbrField = nbrMesh().lookupObject<volScalarField>
        (
            //presume same field name as on this side
            this->dimensionedInternalField().name()
        );

    fieldNbrToOwn = interpolateFromNbrField<scalar>
        (
            nbrPatch().patchField<volScalarField, scalar>(nbrField)
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
    os.writeKeyword("k") << kName_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourRegionName") << neighbourRegionName_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourPatchName") << neighbourPatchName_ 
        << token::END_STATEMENT << nl;        
    os.writeKeyword("relax") << relax_ << token::END_STATEMENT << nl;
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
    fvPatchScalarField,
    interfaceCoupledPressureValue
);

} // End namespace Foam

// ************************************************************************* //
