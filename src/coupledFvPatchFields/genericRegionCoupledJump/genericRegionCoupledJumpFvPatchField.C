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

#include "genericRegionCoupledJumpFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
genericRegionCoupledJumpFvPatchField<Type>::genericRegionCoupledJumpFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    interfaceToInterfaceCoupleManager(p),
    psiName_(this->dimensionedInternalField().name()),
    kName_("k"),
    neighbourRegionName_(),
    neighbourPatchName_(), 
    relax_(1.0),
    nonOrthCorr_(false),
    secondOrder_(false)
{}

template<class Type>
genericRegionCoupledJumpFvPatchField<Type>::genericRegionCoupledJumpFvPatchField
(
    const genericRegionCoupledJumpFvPatchField& grcj,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(grcj, p, iF, mapper),
    interfaceToInterfaceCoupleManager(p),
    psiName_(grcj.psiName_),
    kName_(grcj.kName_),
    neighbourRegionName_(grcj.neighbourRegionName_),
    neighbourPatchName_(grcj.neighbourPatchName_),
    relax_(grcj.relax_),
    nonOrthCorr_(grcj.nonOrthCorr_),
    secondOrder_(grcj.secondOrder_)
{}

template<class Type>
genericRegionCoupledJumpFvPatchField<Type>::genericRegionCoupledJumpFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
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
    relax_(dict.lookupOrDefault<scalar>("relax",1.0)),
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
genericRegionCoupledJumpFvPatchField<Type>::genericRegionCoupledJumpFvPatchField
(
    const genericRegionCoupledJumpFvPatchField& grcj,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(grcj, iF),
    interfaceToInterfaceCoupleManager(grcj),
    psiName_(grcj.psiName_),
    kName_(grcj.kName_),
    neighbourRegionName_(grcj.neighbourRegionName_),
    neighbourPatchName_(grcj.neighbourPatchName_),
    relax_(grcj.relax_),
    nonOrthCorr_(grcj.nonOrthCorr_),
    secondOrder_(grcj.secondOrder_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return patch-normal gradient
template<class Type>
tmp<Field<Type> > genericRegionCoupledJumpFvPatchField<Type>::snGrad() const
{
    notImplemented
    (
        "genericRegionCoupledJumpFvPatchFieldFvPatchFieldTemplates.C\n"
        "tmp<Field<Type> > genericRegionCoupledJumpFvPatchField::snGrad()\n"
        "not implemented"
    );

    return (*this * 0);
}

template<class Type>
tmp<Field<Type> > genericRegionCoupledJumpFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    notImplemented
    (
        "genericRegionCoupledJumpFvPatchFieldFvPatchFieldTemplates.C\n"
        "tmp<Field<Type> > genericRegionCoupledJumpFvPatchField::gradientBoundaryCoeffs()\n"
        "not implemented"
    );

    return (*this * 0);
}

//- Update the patch field coefficients
template<class Type>
void genericRegionCoupledJumpFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Lookup neighbouring patch field
    const GeometricField<Type, fvPatchField, volMesh>& 
        nbrField = nbrMesh().lookupObject<GeometricField<Type, fvPatchField, volMesh>>
        (
            //presume same field name as on this side
            psiName_
        );

    // Calculate interpolated patch field
    Field<Type> fieldNbrToOwn = interpolateFromNbrField<Type>
    (
        nbrPatch().patchField<GeometricField<Type, fvPatchField, volMesh>, Type>(nbrField)
    );

    // Add interfacial pressure jump
    fieldNbrToOwn += valueJump();

    Field<Type>& patchField = *this;

    patchField = *this + relax_*(fieldNbrToOwn - *this);

    updatePhi();

    fixedValueFvPatchField<Type>::updateCoeffs();
}


//- Return the patch flux
template<class Type>
tmp<Field<Type> > genericRegionCoupledJumpFvPatchField<Type>::flux() const
{
    
    dimensionedScalar k
    (
        this->db().time().objectRegistry::lookupObject<IOdictionary>("transportProperties")
        .subDict(refPatch().boundaryMesh().mesh().name()).lookup(kName_)
    );

    return (this->snGrad()*k.value());
}

//- Return the raw coupled patch residual
template<class Type>
scalarField genericRegionCoupledJumpFvPatchField<Type>::rawResidual() const
{
    // Lookup neighbouring patch field
    const GeometricField<Type, fvPatchField, volMesh>& 
        nbrField = nbrMesh().lookupObject<GeometricField<Type, fvPatchField, volMesh>>
        (
            //presume same field name as on this side
            psiName_
        );

    // Calculate interpolated patch field
    Field<Type> fieldNbrToOwn = interpolateFromNbrField<Type>
        (
            nbrPatch().patchField<GeometricField<Type, fvPatchField, volMesh>, Type>(nbrField)
        );
    
    const Field<Type>& fown = *this;

    // Calculate the raw residual
    const scalarField& rawResidual =
        (
            mag(mag(fown - fieldNbrToOwn) - mag(valueJump()))
        );

    return rawResidual;
}

//- Return the normalized coupled patch residual
template<class Type>
scalarField genericRegionCoupledJumpFvPatchField<Type>::normResidual() const
{
    // Lookup neighbouring patch field
    const GeometricField<Type, fvPatchField, volMesh>& 
        nbrField = nbrMesh().lookupObject<GeometricField<Type, fvPatchField, volMesh>>
        (
            //presume same field name as on this side
            psiName_
        );

    // Calculate interpolated patch field
    Field<Type> fieldNbrToOwn = interpolateFromNbrField<Type>
        (
            nbrPatch().patchField<GeometricField<Type, fvPatchField, volMesh>, Type>(nbrField)
        );
    
    const Field<Type>& fown = *this;

    //Calculate normalisation factor
    const scalar n = max
        (
            min
            (
                gMax(mag(fown)), 
                gMax(mag(fieldNbrToOwn + valueJump()))
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
scalar genericRegionCoupledJumpFvPatchField<Type>::ofNormResidual() const
{
    // Lookup neighbouring patch field
    const GeometricField<Type, fvPatchField, volMesh>& 
        nbrField = nbrMesh().lookupObject<GeometricField<Type, fvPatchField, volMesh>>
        (
            //presume same field name as on this side
            psiName_
        );

    // Calculate interpolated patch field
    Field<Type> fieldNbrToOwn = interpolateFromNbrField<Type>
        (
            nbrPatch().patchField<GeometricField<Type, fvPatchField, volMesh>, Type>(nbrField)
        );

    const Field<Type>& fown = *this;

    //Calculate normalisation factor similar to linear system solver
    scalarField jumpRef(refPatch().size(), gAverage(fown) - gAverage(fieldNbrToOwn));

    const scalar n = max
        (
            gSum
            (
                mag((fown-fieldNbrToOwn) - jumpRef)
              + mag(valueJump() - (fown-fieldNbrToOwn))
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
scalar genericRegionCoupledJumpFvPatchField<Type>::maxNormResidual() const
{
    scalar maxNormResidual = gMax(normResidual());

    return maxNormResidual;
}

//- Write
template<class Type>
void genericRegionCoupledJumpFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("k") << kName_ << token::END_STATEMENT << nl;
    os.writeKeyword("relax") << relax_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourRegionName") << neighbourRegionName_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourPatchName") << neighbourPatchName_ 
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
// #   include "genericRegionCoupledJumpFvPatchFieldTemplates.C"
// #endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
