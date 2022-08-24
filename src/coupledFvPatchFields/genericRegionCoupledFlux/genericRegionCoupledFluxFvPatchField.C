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

#include "genericRegionCoupledFluxFvPatchField.H"
#include "genericRegionCoupledJumpFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
genericRegionCoupledFluxFvPatchField<Type>::genericRegionCoupledFluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    interfaceToInterfaceCoupleManager(p),
    psiName_(this->dimensionedInternalField().name()),
    kName_("k"),
    neighbourRegionName_(),
    neighbourPatchName_(),
    nonOrthCorr_(false),
    secondOrder_(false)
{}

template<class Type>
genericRegionCoupledFluxFvPatchField<Type>::genericRegionCoupledFluxFvPatchField
(
    const genericRegionCoupledFluxFvPatchField& grcf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchField<Type>(grcf, p, iF, mapper),
    interfaceToInterfaceCoupleManager(p),
    psiName_(grcf.psiName_),
    kName_(grcf.kName_),
    neighbourRegionName_(grcf.neighbourRegionName_),
    neighbourPatchName_(grcf.neighbourPatchName_),
    nonOrthCorr_(grcf.nonOrthCorr_),
    secondOrder_(grcf.secondOrder_)
{}

template<class Type>
genericRegionCoupledFluxFvPatchField<Type>::genericRegionCoupledFluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    interfaceToInterfaceCoupleManager(p, dict),
    psiName_(this->dimensionedInternalField().name()),
    kName_(dict.lookupOrDefault<word>("k", word::null)),
    neighbourRegionName_
    (
        dict.lookupOrDefault<word>("neighbourRegionName", word::null)
    ),
    neighbourPatchName_
    (
        dict.lookupOrDefault<word>("neighbourPatchName", word::null)
    ),
    nonOrthCorr_(dict.lookupOrDefault<Switch>("nonOrthCorr",false)),
    secondOrder_(dict.lookupOrDefault<Switch>("secondOrder",false))
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        this->evaluate();
    }
}

template<class Type>
genericRegionCoupledFluxFvPatchField<Type>::genericRegionCoupledFluxFvPatchField
(
    const genericRegionCoupledFluxFvPatchField& grcf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(grcf, iF),
    interfaceToInterfaceCoupleManager(grcf),
    psiName_(grcf.psiName_),
    kName_(grcf.kName_),
    neighbourRegionName_(grcf.neighbourRegionName_),
    neighbourPatchName_(grcf.neighbourPatchName_),
    nonOrthCorr_(grcf.nonOrthCorr_),
    secondOrder_(grcf.secondOrder_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void genericRegionCoupledFluxFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    notImplemented
    (
        "genericRegionCoupledFluxFvPatchFieldFvPatchFieldTemplates.C\n"
        "tmp<Field<Type> > genericRegionCoupledFluxFvPatchField::snGrad()\n"
        "not implemented"
    );
}

//- Update the patch field coefficients
template<class Type>
void genericRegionCoupledFluxFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Update the region interface
    const_cast<regionInterface&>(rgInterface()).update();

    // Lookup neighbouring patch field
    const GeometricField<Type, fvPatchField, volMesh>& 
        nbrField = nbrMesh().lookupObject<GeometricField<Type, fvPatchField, volMesh>>
        (
            //presume same field name as on this side
            psiName_
        );

    // Interpolate flux face values from neighbour patch
    tmp<Field<Type>> tnbrFlux = 
        refCast<const genericRegionCoupledJumpFvPatchField<Type>>
        (
            nbrPatch()
            .patchField<GeometricField<Type, fvPatchField, volMesh>, Type>(nbrField)
        ).flux();

    const Field<Type>& nbrFlux = tnbrFlux();

    // Calculate interpolated patch field
    Field<Type> fluxNbrToOwn = interpolateFromNbrField<Type>(nbrFlux);

    // Enforce flux matching
    fluxNbrToOwn *= -1.0; 
    fluxNbrToOwn += fluxJump();

    dimensionedScalar k
    (
        this->db().time().objectRegistry::lookupObject<IOdictionary>("transportProperties")
        .subDict(refPatch().boundaryMesh().mesh().name()).lookup(kName_)
    );

    // Add interfacial pressure Flux
    this->gradient() = fluxNbrToOwn/k.value();

    fixedGradientFvPatchField<Type>::updateCoeffs();
}

//- Return the raw coupled patch residual
template<class Type>
scalarField genericRegionCoupledFluxFvPatchField<Type>::rawResidual() const
{
    // Lookup neighbouring patch field
    const GeometricField<Type, fvPatchField, volMesh>& 
        nbrField = nbrMesh().lookupObject<GeometricField<Type, fvPatchField, volMesh>>
        (
            //presume same field name as on this side
            psiName_
        );

    // Interpolate flux face values from neighbour patch
    tmp<Field<Type>> tnbrFlux = 
        refCast<const genericRegionCoupledJumpFvPatchField<Type>>
        (
            nbrPatch()
            .patchField<GeometricField<Type, fvPatchField, volMesh>, Type>(nbrField)
        ).flux();

    const Field<Type>& nbrFlux = tnbrFlux();

    // Calculate interpolated patch field
    Field<Type> fluxNbrToOwn = interpolateFromNbrField<Type>(nbrFlux);

    dimensionedScalar k
    (
        this->db().time().objectRegistry::lookupObject<IOdictionary>("transportProperties")
        .subDict(refPatch().boundaryMesh().mesh().name()).lookup(kName_)
    );

    const Field<Type> fluxOwn = this->snGrad()*k.value();

    const tmp<Field<Type>> tmpFluxJump = fluxJump();
    const Field<Type>& fluxJump =  tmpFluxJump();

    Field<Type> residual = (fluxOwn + fluxNbrToOwn) - fluxJump;

    // Calculate the raw residual
    const tmp<scalarField> tmpRawResidual = mag(residual);
    const scalarField rawResidual = tmpRawResidual();
    
    // Info<< nl
    // << psiName_ << " fluxOwn:" << nl
    // << "  max: " << gMax(fluxOwn) << nl
    // << "  min: " << gMin(fluxOwn) << nl
    // << "  mean: " << gAverage(fluxOwn) << nl
    // << psiName_ << " fluxNbrToOwn:" << nl
    // << "  max: " << gMax(fluxNbrToOwn) << nl
    // << "  min: "<< gMax(fluxNbrToOwn) << nl
    // << "  mean: " << gAverage(fluxNbrToOwn) << nl
    // << psiName_ << " fluxJump:" << nl
    // << "  max: " << gMax(fluxJump) << nl
    // << "  min: " << gMin(fluxJump) << nl
    // << "  mean: " << gAverage(fluxJump) << nl
    // << psiName_ << " residual:" << nl
    // << "  max: " << gMax(residual) << nl
    // << "  min: " << gMin(residual) << nl
    // << "  mean: " << gAverage(residual) << nl
    // << psiName_ << " rawResidual:" << nl
    // << "  max: " << gMax(rawResidual) << nl
    // << "  min: " << gMin(rawResidual) << nl
    // << "  mean: " << gAverage(rawResidual) << nl
    // << endl;

    return rawResidual;
}

//- Return the normalized coupled patch residual
template<class Type>
scalarField genericRegionCoupledFluxFvPatchField<Type>::normResidual() const
{
    // Lookup neighbouring patch field
    const GeometricField<Type, fvPatchField, volMesh>& 
        nbrField = nbrMesh().lookupObject<GeometricField<Type, fvPatchField, volMesh>>
        (
            //presume same field name as on this side
            psiName_
        );

    // Interpolate flux face values from neighbour patch
    tmp<Field<Type>> tnbrFlux = 
        refCast<const genericRegionCoupledJumpFvPatchField<Type>>
        (
            nbrPatch()
            .patchField<GeometricField<Type, fvPatchField, volMesh>, Type>(nbrField)
        ).flux();

    const Field<Type>& nbrFlux = tnbrFlux();

    // Calculate interpolated patch field
    Field<Type> fluxNbrToOwn = interpolateFromNbrField<Type>(nbrFlux);

    dimensionedScalar k
    (
        this->db().time().objectRegistry::lookupObject<IOdictionary>("transportProperties")
        .subDict(refPatch().boundaryMesh().mesh().name()).lookup(kName_)
    );

    const Field<Type> fluxOwn = this->snGrad()*k.value();

    //Calculate normalisation factor
    const scalar n = max
        (
            min
            (
                gMax(mag(fluxOwn)), 
                //gMax(mag(fluxNbrToOwn + fluxJump()))
                gMax(mag(fluxNbrToOwn))
            ),
            SMALL
        );

    //Return normalised residual
    return 
        (
            rawResidual()/n
        );
}

//- Return the normalized coupled patch residual
//- normalised similar to linear system solver residuals
template<class Type>
scalar genericRegionCoupledFluxFvPatchField<Type>::ofNormResidual() const
{
    // Lookup neighbouring patch field
    const GeometricField<Type, fvPatchField, volMesh>& 
        nbrField = nbrMesh().lookupObject<GeometricField<Type, fvPatchField, volMesh>>
        (
            //presume same field name as on this side
            psiName_
        );

    // Interpolate flux face values from neighbour patch
    tmp<Field<Type>> tnbrFlux = 
        refCast<const genericRegionCoupledJumpFvPatchField<Type>>
        (
            nbrPatch()
            .patchField<GeometricField<Type, fvPatchField, volMesh>, Type>(nbrField)
        ).flux();

    const Field<Type>& nbrFlux = tnbrFlux();

    // Calculate interpolated patch field
    Field<Type> fluxNbrToOwn = interpolateFromNbrField<Type>(nbrFlux);

    dimensionedScalar k
    (
        this->db().time().objectRegistry::lookupObject<IOdictionary>("transportProperties")
        .subDict(refPatch().boundaryMesh().mesh().name()).lookup(kName_)
    );

    const Field<Type> fluxOwn = this->snGrad()*k.value();

    //Calculate normalisation factor similar to linear system solver
    const Field<Type> fluxRef(refPatch().size(), gAverage(fluxOwn) - gAverage(fluxNbrToOwn));

    const scalar n = max
        (
            gSum
            (
                mag((fluxOwn-fluxNbrToOwn) - fluxRef)
              + mag(fluxJump() - (fluxOwn-fluxNbrToOwn))
            ),
            SMALL
        );

    return 
        (
            gSum(rawResidual())/n
        );

}

//- Return the maximum normalized coupled patch residual
template<class Type>
scalar genericRegionCoupledFluxFvPatchField<Type>::maxNormResidual() const
{
    scalar maxNormResidual = gMax(normResidual());

    return maxNormResidual;
}

//- Write
template<class Type>
void genericRegionCoupledFluxFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("k") << kName_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourRegionName") << neighbourRegionName_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourPatchName") << neighbourPatchName_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourFieldName") << psiName_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("nonOrthCorr") << nonOrthCorr_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("secondOrder") << secondOrder_ 
        << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "genericRegionCoupledFluxFvPatchFieldTemplates.C"
// #endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
