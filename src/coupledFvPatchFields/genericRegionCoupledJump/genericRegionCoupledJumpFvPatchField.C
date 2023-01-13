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
    neighbourRegionName_(),
    neighbourPatchName_(),
    neighbourFieldName_(),
    kName_("k"),
    KName_("K"),
    relaxModel_(p.boundaryMesh().mesh().time()),
    nonOrthCorr_(false),
    secondOrder_(false)
{
    relaxModel_.initialize(*this);
}

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
    neighbourRegionName_(grcj.neighbourRegionName_),
    neighbourPatchName_(grcj.neighbourPatchName_),
    neighbourFieldName_(grcj.neighbourFieldName_),
    kName_(grcj.kName_),
    KName_(grcj.KName_),
    relaxModel_(grcj.relaxModel_),
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
    neighbourRegionName_(dict.lookup("neighbourRegionName")),
    neighbourPatchName_(dict.lookup("neighbourPatchName")),
    neighbourFieldName_(this->dimensionedInternalField().name()),
    kName_(dict.lookup("k")),
    KName_(dict.lookupOrDefault<word>("K", word::null)),
    relaxModel_(p.boundaryMesh().mesh().time(), dict),
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

    relaxModel_.initialize(*this);
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
    neighbourRegionName_(grcj.neighbourRegionName_),
    neighbourPatchName_(grcj.neighbourPatchName_),
    neighbourFieldName_(grcj.neighbourFieldName_),
    kName_(grcj.kName_),
    KName_(grcj.KName_),
    relaxModel_(grcj.relaxModel_),
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
tmp<Field<Type> > 
genericRegionCoupledJumpFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    notImplemented
    (
        "genericRegionCoupledJumpFvPatchFieldFvPatchFieldTemplates.C\n"
        "tmp<Field<Type> > genericRegionCoupledJumpFvPatchField::gradientBoundaryCoeffs()\n"
        "not implemented"
    );

    return (*this * 0);
}

template<class Type>
void genericRegionCoupledJumpFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Update and correct the region interface physics
    const_cast<regionInterface&>(rgInterface()).update();
    const_cast<regionInterface&>(rgInterface()).correct();

    // Lookup neighbouring patch field
    const GeometricField<Type, fvPatchField, volMesh>& nbrField =
        nbrMesh().lookupObject<GeometricField<Type, fvPatchField, volMesh> >
        (
            // same field name as on this side
            this->dimensionedInternalField().name()
        );

    // Calculate interpolated patch field
    Field<Type> fieldNbrToOwn = interpolateFromNbrField<Type>
    (
        nbrPatch()
        .patchField<GeometricField<Type, fvPatchField, volMesh>, Type>(nbrField)
    );

    // Add interfacial jump
    fieldNbrToOwn += valueJump();

    // Enforce fixed value condition
    fvPatchField<Type>::operator==
    (
        fieldNbrToOwn
    );

    // Relax fixed value condition
    relaxModel_.relax(*this);

    updatePhi();

    fixedValueFvPatchField<Type>::updateCoeffs();
}


//- Return the patch flux
template<class Type>
tmp<Field<Type> > genericRegionCoupledJumpFvPatchField<Type>::flux() const
{
    // Get the diffusivity
    scalarField k(this->patch().size(), pTraits<scalar>::zero);

    if ( this->db().objectRegistry::foundObject<volScalarField>(kName_) )
    {
        k = this->patch().template lookupPatchField<volScalarField, scalar>(kName_);
    }
    else
    {
        k = dimensionedScalar
        (
            this->db().objectRegistry::
            lookupObject<IOdictionary>("transportProperties")
            .lookup(kName_)
        ).value();
    }

    return (this->snGrad()*k);
}

template<class Type>
scalarField genericRegionCoupledJumpFvPatchField<Type>::rawResidual() const
{
    // Lookup neighbouring patch field
    const GeometricField<Type, fvPatchField, volMesh>& nbrField =
        nbrMesh().lookupObject<GeometricField<Type, fvPatchField, volMesh> >
        (
            //same field name as on this side
            this->dimensionedInternalField().name()
        );

    // Calculate interpolated patch field
    Field<Type> fieldNbrToOwn =
        interpolateFromNbrField<Type>
        (
            nbrPatch()
            .patchField<GeometricField<Type, fvPatchField, volMesh>, Type>(nbrField)
        );
    
    const Field<Type>& fown = *this;

    const tmp<Field<Type>> tmpValueJump = valueJump();
    const Field<Type>& valueJump =  tmpValueJump();

    Field<Type> residual = (fown - fieldNbrToOwn) - valueJump;

    // Calculate the raw residual
    const tmp<scalarField> tmpRawResidual = mag(residual);
    const scalarField rawResidual = tmpRawResidual();

    if (debug)
    {
        Info<< nl
        << this->dimensionedInternalField().name() << " fown:" << nl
        << "  max: " << gMax(fown) << nl
        << "  min: " << gMin(fown) << nl
        << "  mean: " << gAverage(fown) << nl
        << this->dimensionedInternalField().name() << " fieldNbrToOwn:" << nl
        << "  max: " << gMax(fieldNbrToOwn) << nl
        << "  min: "<< gMax(fieldNbrToOwn) << nl
        << "  mean: " << gAverage(fieldNbrToOwn) << nl
        << this->dimensionedInternalField().name() << " valueJump:" << nl
        << "  max: " << gMax(valueJump) << nl
        << "  min: " << gMin(valueJump) << nl
        << "  mean: " << gAverage(valueJump)<< nl
        << this->dimensionedInternalField().name() << " residual:" << nl
        << "  max: " << gMax(residual) << nl
        << "  min: " << gMin(residual) << nl
        << "  mean: " << gAverage(residual) << endl
        << this->dimensionedInternalField().name() << " rawResidual:" << nl
        << "  max: " << gMax(rawResidual) << nl
        << "  min: " << gMin(rawResidual) << nl
        << "  mean: " << gAverage(rawResidual) << nl
        << endl;
    }

    return rawResidual;
}

template<class Type>
scalarField genericRegionCoupledJumpFvPatchField<Type>::normResidual() const
{
    // Lookup neighbouring patch field
    const GeometricField<Type, fvPatchField, volMesh>& nbrField =
        nbrMesh().lookupObject<GeometricField<Type, fvPatchField, volMesh>>
        (
            // same field name as on this side
            this->dimensionedInternalField().name()
        );

    // Calculate interpolated patch field
    Field<Type> fieldNbrToOwn = interpolateFromNbrField<Type>
    (
        nbrPatch()
        .patchField<GeometricField<Type, fvPatchField, volMesh>, Type>(nbrField)
    );
    
    const Field<Type>& fown = *this;

    const tmp<Field<Type>> tmpValueJump = valueJump();
    const Field<Type>& valueJump =  tmpValueJump();

    //Calculate normalisation factor
    const scalar n =
        max
        (
            min
            (
                gMax(mag(fown)), 
                gMax(mag(fieldNbrToOwn + valueJump))
            ),
            SMALL
        );

    //Return normalised residual
    return 
    (
        rawResidual()/n
    );
}

template<class Type>
scalar genericRegionCoupledJumpFvPatchField<Type>::ofNormResidual() const
{
    // Lookup neighbouring patch field
    const GeometricField<Type, fvPatchField, volMesh>& nbrField =
        nbrMesh().lookupObject<GeometricField<Type, fvPatchField, volMesh> >
        (
            //same field name as on this side
            this->dimensionedInternalField().name()
        );

    // Calculate interpolated patch field
    Field<Type> fieldNbrToOwn = interpolateFromNbrField<Type>
    (
        nbrPatch()
        .patchField<GeometricField<Type, fvPatchField, volMesh>, Type>(nbrField)
    );

    const Field<Type>& fown = *this;

    //Calculate normalisation factor similar to linear system solver
    const Field<Type> jumpRef
    (
        refPatch().size(),
        gAverage(fown) - gAverage(fieldNbrToOwn)
    );

    const scalar n =
        max
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

template<class Type>
scalar genericRegionCoupledJumpFvPatchField<Type>::maxNormResidual() const
{
    scalar maxNormResidual = gMax(normResidual());

    return maxNormResidual;
}

template<class Type>
void genericRegionCoupledJumpFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("neighbourRegionName") << neighbourRegionName_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourPatchName") << neighbourPatchName_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourFieldName") << neighbourFieldName_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("k") << kName_ << token::END_STATEMENT << nl;
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
