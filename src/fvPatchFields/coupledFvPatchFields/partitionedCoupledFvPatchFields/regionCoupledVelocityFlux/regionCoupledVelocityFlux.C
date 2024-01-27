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

#include "regionCoupledVelocityFlux.H"
#include "regionCoupledVelocityValue.H"
#include "addToRunTimeSelectionTable.H"
#include "primitiveFieldsFwd.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupledVelocityFlux::
regionCoupledVelocityFlux
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    genericRegionCoupledFluxFvPatchField<vector>(p, iF)
{}


Foam::regionCoupledVelocityFlux::
regionCoupledVelocityFlux
(
    const regionCoupledVelocityFlux& icvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    genericRegionCoupledFluxFvPatchField<vector>(icvf, p, iF, mapper)
{}


Foam::regionCoupledVelocityFlux::
regionCoupledVelocityFlux
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    genericRegionCoupledFluxFvPatchField<vector>(p, iF, dict)
{}


Foam::regionCoupledVelocityFlux::
regionCoupledVelocityFlux
(
    const regionCoupledVelocityFlux& icvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    genericRegionCoupledFluxFvPatchField<vector>(icvf, iF)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> regionCoupledVelocityFlux::fluxJump() const
{
    const vectorField& nf = capInterface().aMesh().faceAreaNormals();

    // Lookup neighbouring patch field
    const volVectorField& nbrUField =
        nbrMesh().lookupObject<volVectorField>
        (
            // same field name as on this side
            this->dimensionedInternalField().name()
        );

    // Get flux face values from neighbour patch
    tmp<vectorField> tnbrUFlux =
        refCast<const genericRegionCoupledJumpFvPatchField<vector>>
        (
            nbrPatch()
            .patchField<volVectorField, vector>(nbrUField)
        ).flux();

    const vectorField& nbrUFlux = tnbrUFlux();

    // Calculate interpolated patch field
    vectorField UfluxNbrToOwn = interpolateFromNbrField<vector>(nbrUFlux);

    // Enforce flux matching
    UfluxNbrToOwn *= -1.0;

    // surface tension
    const areaScalarField& sigma = capInterface().sigma();

    // surface velocity terms
    const areaVectorField& Us = capInterface().Us();

    areaScalarField divSU = fac::div(Us);
    divSU.correctBoundaryConditions();

    areaTensorField gradSU = fac::grad(Us);

    vectorField surfaceTensionForce =
        sigma
       *fac::edgeIntegrate
        (
            capInterface().aMesh().Le()
            *capInterface().aMesh().edgeLengthCorrection()
        )().internalField();

    vectorField tangentialSurfaceTensionForce =
        surfaceTensionForce
      - sigma
       *capInterface().aMesh().faceCurvatures().internalField()*nf;

    dimensionedScalar muFluidNbr
    (
        nbrMesh().lookupObject<IOdictionary>("transportProperties")
        .lookup("mu")
    );

    dimensionedScalar muFluid
    (
        refMesh().lookupObject<IOdictionary>("transportProperties")
        .lookup("mu")
    );

    return
    (
      - nf*(nf & UfluxNbrToOwn)
      + tangentialSurfaceTensionForce
      - muFluid.value()*nf*divSU.internalField()
      + (muFluidNbr.value() - muFluid.value())
       *(gradSU.internalField()&nf)
    );
}

const regionInterfaces::capillaryInterface&
regionCoupledVelocityFlux::capInterface() const
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
    fvPatchVectorField,
    regionCoupledVelocityFlux
);

} // End namespace Foam

// ************************************************************************* //
