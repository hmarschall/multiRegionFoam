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

#include "primitiveFieldsFwd.H"
#include "regionCoupledPotentialJump.H"
#include "addToRunTimeSelectionTable.H"
#include "currentTransferInterface.H"

#include "patchToPatchInterpolation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupledPotentialJump::
regionCoupledPotentialJump
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    genericRegionCoupledJumpFvPatchField<scalar>(p, iF)
{}


Foam::regionCoupledPotentialJump::
regionCoupledPotentialJump
(
    const regionCoupledPotentialJump& icpv,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    genericRegionCoupledJumpFvPatchField<scalar>(icpv, p, iF, mapper)
{}


Foam::regionCoupledPotentialJump::
regionCoupledPotentialJump
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    genericRegionCoupledJumpFvPatchField<scalar>(p, iF, dict)
{}


Foam::regionCoupledPotentialJump::
regionCoupledPotentialJump
(
    const regionCoupledPotentialJump& icpv,
    const DimensionedField<scalar, volMesh>& iF
)
:
    genericRegionCoupledJumpFvPatchField<scalar>(icpv, iF)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<scalarField> regionCoupledPotentialJump::valueJump() const
{
    // Lookup neighbouring patch field
    const GeometricField<scalar, fvPatchField, volMesh>& nbrField =
        nbrMesh().lookupObject<GeometricField<scalar, fvPatchField, volMesh>>
        (
            // same field name as on this side
            this->dimensionedInternalField().name()
        );

    // Calculate interpolated patch field
    Field<scalar> fieldNbrToOwn = interpolateFromNbrField<scalar>
    (
        nbrPatch()
        .patchField<GeometricField<scalar, fvPatchField, volMesh>, scalar>(nbrField)
    );

    // Get the jump field from the interface
    const word fluidPhaseName = currentTransInterface().fluidPhaseName();
    const word fluidPatchName = currentTransInterface().fluidPatchName();
    const scalar etaSign = currentTransInterface().etaSign();

    // Get the values of the potential jump
    const regionType& fluidPhase = this->db().time().objectRegistry::lookupObject<regionType>(currentTransInterface().fluidPhaseName() + "Dict");
    const volScalarField& jump = fluidPhase.mesh().lookupObject<volScalarField>("jump");
    label fluidPatchID = jump.mesh().boundaryMesh().findPatchID(fluidPatchName);
    const fvPatchField<scalar> jumpPatch = jump.boundaryField()[fluidPatchID];

    // Map it onto the current boundary of this mesh
    const polyPatch& fluidPatch = jump.mesh().boundaryMesh()[fluidPatchID];
    patchToPatchInterpolation interpolator2(fluidPatch, patch().patch());  // ownPatch);
    scalarField jumpFluidToOwn = interpolator2.faceInterpolate(jumpPatch);

    return (etaSign*jumpFluidToOwn);
}

const regionInterfaces::currentTransferInterface&
regionCoupledPotentialJump::currentTransInterface() const
{
    if(   rgInterface().type()
       != regionInterfaces::currentTransferInterface::typeName )
    {
        FatalErrorInFunction
            << this->typeName << " BC can only "
            << "be used in combination with a "
            << regionInterfaces::currentTransferInterface::typeName
            << endl
            << exit(FatalError);
    }

    return refCast<const regionInterfaces::currentTransferInterface>
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
    regionCoupledPotentialJump
);

} // End namespace Foam

// ************************************************************************* //
