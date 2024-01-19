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

#include "genericInterfaceCoupledPressureFlux.H"
#include "genericInterfaceCoupledPressureValue.H"
#include "addToRunTimeSelectionTable.H"
#include "primitiveFieldsFwd.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::genericInterfaceCoupledPressureFlux::
genericInterfaceCoupledPressureFlux
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    genericRegionCoupledFluxFvPatchField<scalar>(p, iF)
{}


Foam::genericInterfaceCoupledPressureFlux::
genericInterfaceCoupledPressureFlux
(
    const genericInterfaceCoupledPressureFlux& icpf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    genericRegionCoupledFluxFvPatchField<scalar>(icpf, p, iF, mapper)
{}


Foam::genericInterfaceCoupledPressureFlux::
genericInterfaceCoupledPressureFlux
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    genericRegionCoupledFluxFvPatchField<scalar>(p, iF, dict)
{}


Foam::genericInterfaceCoupledPressureFlux::
genericInterfaceCoupledPressureFlux
(
    const genericInterfaceCoupledPressureFlux& icpf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    genericRegionCoupledFluxFvPatchField<scalar>(icpf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> genericInterfaceCoupledPressureFlux::fluxJump() const
{
    // const volVectorField& U =
    //     refMesh().objectRegistry::lookupObject<volVectorField>("U");

    // const volVectorField& nbrU =
    //     nbrMesh().lookupObject<volVectorField>("U");

    // const surfaceScalarField& phi =
    //     refMesh().objectRegistry::lookupObject<surfaceScalarField>("phi");

    // const surfaceScalarField& nbrPhi =
    //     nbrMesh().lookupObject<surfaceScalarField>("phi");

    // vectorField nbrLaplacianU = interpolateFromNbrField<vector>
    //     (
    //         fvc::laplacian(nbrU)()
    //         .boundaryField()[nbrPatch().index()]
    //         .patchInternalField()
    //     );

    // vectorField nbrDDtU = interpolateFromNbrField<vector>
    //     (
    //         fvc::DDt(nbrPhi, nbrU)()
    //         .boundaryField()[nbrPatch().index()]
    //         .patchInternalField()
    //     );

    // dimensionedScalar muFluidNbr
    // (
    //     nbrMesh().lookupObject<IOdictionary>("transportProperties")
    //     .lookup("mu")
    // );

    // dimensionedScalar muFluid
    // (
    //     refMesh().lookupObject<IOdictionary>("transportProperties")
    //     .lookup("mu")
    // );

    // dimensionedScalar rhoFluidNbr
    // (
    //     nbrMesh().lookupObject<IOdictionary>("transportProperties")
    //     .lookup("rho")
    // );

    // dimensionedScalar rhoFluid
    // (
    //     refMesh().lookupObject<IOdictionary>("transportProperties")
    //     .lookup("rho")
    // );

    // return
    // (
    //     (
    //         muFluid.value()/rhoFluid.value()
    //         *(
    //             fvc::laplacian(U)()
    //             .boundaryField()[this->patch().index()]
    //             .patchInternalField()
    //         )
    //       - muFluidNbr.value()/rhoFluidNbr.value()*
    //         nbrLaplacianU
    //     )
    //   + (
    //         nbrDDtU
    //       - (
    //             fvc::DDt(phi,U)()
    //             .boundaryField()[this->patch().index()]
    //             .patchInternalField()
    //         )
    //     )
    // ) & patch().nf();

    // Lookup neighbouring patch field
    const volScalarField& nbrPField =
        nbrMesh().lookupObject<volScalarField>
        (
            // same field name as on this side
            this->dimensionedInternalField().name()
        );

    // Interpolate flux face values from neighbour patch
    tmp<scalarField> tnbrPFlux =
        refCast<const genericRegionCoupledJumpFvPatchField<scalar>>
        (
            nbrPatch()
            .patchField<volScalarField, scalar>(nbrPField)
        ).flux();

    const scalarField& nbrPFlux = tnbrPFlux();

    // Calculate interpolated patch field
    scalarField pfluxNbrToOwn = interpolateFromNbrField<scalar>(nbrPFlux);

    // Enforce flux matching
    pfluxNbrToOwn *= -1.0;

    dimensionedScalar rhoFluid
    (
        refMesh().lookupObject<IOdictionary>("transportProperties")
        .lookup("rho")
    );

    vectorField nB = refMesh().boundary()[refPatchID()].nf();

    const volVectorField& U =
        refMesh().objectRegistry::lookupObject<volVectorField>("U");

    return
    (
      - pfluxNbrToOwn
      - rhoFluid.value()
       *(
            nB&fvc::ddt(U)().boundaryField()[refPatchID()]
        )
    );
}

const regionInterfaces::capillaryInterface&
genericInterfaceCoupledPressureFlux::capInterface() const
{
    if(   rgInterface().type()
       != regionInterfaces::capillaryInterface::typeName )
    {
        FatalErrorInFunction
            << this->typeName << " BC can only "
            << "be used in combination with a "
            << regionInterfaces::capillaryInterface::typeName
            << endl
            << exit(FatalError);
    }

    return refCast<const regionInterfaces::capillaryInterface>
        (
            rgInterface()
        );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    genericInterfaceCoupledPressureFlux
);

} // End namespace Foam

// ************************************************************************* //
