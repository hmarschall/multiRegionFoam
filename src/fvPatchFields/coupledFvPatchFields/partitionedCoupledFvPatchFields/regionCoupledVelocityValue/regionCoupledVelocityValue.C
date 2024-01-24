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

#include "regionCoupledVelocityValue.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupledVelocityValue::
regionCoupledVelocityValue
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    genericRegionCoupledJumpFvPatchField<vector>(p, iF)
{}


Foam::regionCoupledVelocityValue::
regionCoupledVelocityValue
(
    const regionCoupledVelocityValue& icvv,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    genericRegionCoupledJumpFvPatchField<vector>(icvv, p, iF, mapper)
{}


Foam::regionCoupledVelocityValue::
regionCoupledVelocityValue
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    genericRegionCoupledJumpFvPatchField<vector>(p, iF, dict)
{}


Foam::regionCoupledVelocityValue::
regionCoupledVelocityValue
(
    const regionCoupledVelocityValue& icvv,
    const DimensionedField<vector, volMesh>& iF
)
:
    genericRegionCoupledJumpFvPatchField<vector>(icvv, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionCoupledVelocityValue::updatePhi()
{
    //- Non-const access to flux on patch
    fvsPatchField<scalar>& patchPhiField = const_cast<fvsPatchField<scalar>& >
    (
        this->db().lookupObject<surfaceScalarField>("phi")
        .boundaryField()[this->patch().index()]
    );

    //- Get the flux on neighboring patch
    surfaceScalarField nbrPhi =
        nbrMesh().lookupObject<surfaceScalarField>("phi");

    //- Impose interpolated flux field
    patchPhiField = interpolateFromNbrField<scalar>
        (
            nbrPatch().patchField<surfaceScalarField, scalar>(nbrPhi)
        )*(-1.); // consider outer normals pointing in opposite directions
}


//- Zero velocity jump
tmp<vectorField> Foam::regionCoupledVelocityValue::valueJump() const
{
    const vectorField nf = refMesh().boundary()[refPatchID()].nf();

    vectorField UsNbrToOwn = interpolateFromNbrField<vector>(capInterface().Us());

    const volVectorField& U =
        refMesh().objectRegistry::lookupObject<volVectorField>("U");

    const volScalarField& rho =
        refMesh().objectRegistry::lookupObject<volScalarField>("rho");

    return
    (
        - nf*(nf & UsNbrToOwn)
        + nf * fvc::meshPhi(rho,U)().boundaryField()[refPatchID()]/
         refMesh().boundary()[refPatchID()].magSf()
    );
}

const regionInterfaces::capillaryInterface&
regionCoupledVelocityValue::capInterface() const
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
    regionCoupledVelocityValue
);

} // End namespace Foam

// ************************************************************************* //

