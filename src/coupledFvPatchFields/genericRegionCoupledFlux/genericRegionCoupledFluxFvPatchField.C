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
    neighbourRegionName_(),
    neighbourPatchName_(),
    neighbourFieldName_(),
    kName_("k"),
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
    neighbourRegionName_(grcf.neighbourRegionName_),
    neighbourPatchName_(grcf.neighbourPatchName_),
    neighbourFieldName_(grcf.neighbourFieldName_),
    kName_(grcf.kName_),
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
    neighbourRegionName_(dict.lookup("neighbourRegionName")),
    neighbourPatchName_(dict.lookup("neighbourPatchName")),
    neighbourFieldName_(this->dimensionedInternalField().name()),
    kName_(dict.lookup("k")),
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

    // Coupled fields should have same names,
    // e.g. there is only one temperature -> T
    // (disambiguous since fields are registered to different meshes)
    if (dict.found("neighbourPatchName"))
    {
        if
        (
            !this->db().objectRegistry::
            foundObject<GeometricField<Type, fvPatchField, volMesh> >
            (dict.lookup("neighbourFieldName"))
        )
        {
            FatalError
                << "\nIncorrect neigbour field name " 
                << dict.lookup("neighbourFieldName")
                << " instead " << this->dimensionedInternalField().name()
                << exit(FatalError);
        }
    }
    else
    {
        FatalError
            << "\nNeigbour field name not found but needed for coupling manager"
            << " Provide neighbourFieldName: " 
            << this->dimensionedInternalField().name()
            << exit(FatalError);
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
    neighbourRegionName_(grcf.neighbourRegionName_),
    neighbourPatchName_(grcf.neighbourPatchName_),
    neighbourFieldName_(grcf.neighbourFieldName_),
    kName_(grcf.kName_),
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

template<class Type>
void genericRegionCoupledFluxFvPatchField<Type>::updateCoeffs()
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

    // Interpolate flux face values from neighbour patch
    tmp<Field<Type> > tnbrFlux = 
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
            refPatch().boundaryMesh().mesh().objectRegistry::
            lookupObject<IOdictionary>("transportProperties")
            .lookup(kName_)
        ).value();
    }

    // Add interfacial flux
    this->gradient() = fluxNbrToOwn/k;

    fixedGradientFvPatchField<Type>::updateCoeffs();
}

template<class Type>
scalarField genericRegionCoupledFluxFvPatchField<Type>::rawResidual() const
{
    // Lookup neighbouring patch field
    const GeometricField<Type, fvPatchField, volMesh>& 
        nbrField = nbrMesh().lookupObject<GeometricField<Type, fvPatchField, volMesh> >
        (
            // same field name as on this side
            this->dimensionedInternalField().name()
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
            refPatch().boundaryMesh().mesh().objectRegistry::
            lookupObject<IOdictionary>("transportProperties")
            .lookup(kName_)
        ).value();
    }

    const Field<Type> fluxOwn = this->snGrad()*k;

    const tmp<Field<Type>> tmpFluxJump = fluxJump();
    const Field<Type>& fluxJump =  tmpFluxJump();

    Field<Type> residual = (fluxOwn + fluxNbrToOwn) - fluxJump;

    // Calculate the raw residual
    const tmp<scalarField> tmpRawResidual = mag(residual);
    const scalarField rawResidual = tmpRawResidual();

    if (debug)
    {
        Info<< nl
        << this->dimensionedInternalField().name() << " fluxOwn:" << nl
        << "  max: " << gMax(fluxOwn) << nl
        << "  min: " << gMin(fluxOwn) << nl
        << "  mean: " << gAverage(fluxOwn) << nl
        << this->dimensionedInternalField().name() << " fluxNbrToOwn:" << nl
        << "  max: " << gMax(fluxNbrToOwn) << nl
        << "  min: "<< gMax(fluxNbrToOwn) << nl
        << "  mean: " << gAverage(fluxNbrToOwn) << nl
        << this->dimensionedInternalField().name() << " fluxJump:" << nl
        << "  max: " << gMax(fluxJump) << nl
        << "  min: " << gMin(fluxJump) << nl
        << "  mean: " << gAverage(fluxJump) << nl
        << this->dimensionedInternalField().name() << " residual:" << nl
        << "  max: " << gMax(residual) << nl
        << "  min: " << gMin(residual) << nl
        << "  mean: " << gAverage(residual) << nl
        << this->dimensionedInternalField().name() << " rawResidual:" << nl
        << "  max: " << gMax(rawResidual) << nl
        << "  min: " << gMin(rawResidual) << nl
        << "  mean: " << gAverage(rawResidual) << nl
        << endl;
    }

    return rawResidual;
}

template<class Type>
scalar genericRegionCoupledFluxFvPatchField<Type>::normResidual() const
{
    // Lookup neighbouring patch field
    const GeometricField<Type, fvPatchField, volMesh>& 
        nbrField = nbrMesh().lookupObject<GeometricField<Type, fvPatchField, volMesh>>
        (
            //same field name as on this side
            this->dimensionedInternalField().name()
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
            refPatch().boundaryMesh().mesh().objectRegistry::
            lookupObject<IOdictionary>("transportProperties")
            .lookup(kName_)
        ).value();
    }

    const Field<Type> fluxOwn = this->snGrad()*k;

    //Calculate normalisation factor
    const scalar n =
        max
        (
            min
            (
                Foam::sqrt(gSum(magSqr(fluxOwn))), 
                Foam::sqrt(gSum(magSqr(-1.0*fluxNbrToOwn + fluxJump())))
                //gMax(mag(fluxNbrToOwn))
            ),
            SMALL
        );

    //Return normalised residual
    return 
    (
        Foam::sqrt(gSum(magSqr(rawResidual())))/n
    );
}

template<class Type>
scalar genericRegionCoupledFluxFvPatchField<Type>::ofNormResidual() const
{
    // Lookup neighbouring patch field
    const GeometricField<Type, fvPatchField, volMesh>& 
        nbrField =
        nbrMesh().lookupObject<GeometricField<Type, fvPatchField, volMesh> >
        (
            // same field name as on this side
            this->dimensionedInternalField().name()
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
            refPatch().boundaryMesh().mesh().objectRegistry::
            lookupObject<IOdictionary>("transportProperties")
            .lookup(kName_)
        ).value();
    }

    const Field<Type> fluxOwn = this->snGrad()*k;

    //Calculate normalisation factor similar to linear system solver
    const Field<Type> fluxRef
    (
        refPatch().size(),
        gAverage(fluxOwn) - gAverage(fluxNbrToOwn)
    );

    const scalar n =
        max
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

template<class Type>
void genericRegionCoupledFluxFvPatchField<Type>::write
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
// #   include "genericRegionCoupledFluxFvPatchFieldTemplates.C"
// #endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
