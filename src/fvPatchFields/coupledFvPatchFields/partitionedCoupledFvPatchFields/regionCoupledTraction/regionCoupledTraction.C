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

#include "regionCoupledTraction.H"
#include "addToRunTimeSelectionTable.H"
#include "primitiveFieldsFwd.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupledTraction::
regionCoupledTraction
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    genericRegionCoupledFluxFvPatchField<vector>(p, iF)
{}


Foam::regionCoupledTraction::
regionCoupledTraction
(
    const regionCoupledTraction& icvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    genericRegionCoupledFluxFvPatchField<vector>(icvf, p, iF, mapper)
{}


Foam::regionCoupledTraction::
regionCoupledTraction
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    genericRegionCoupledFluxFvPatchField<vector>(p, iF, dict)
{}


Foam::regionCoupledTraction::
regionCoupledTraction
(
    const regionCoupledTraction& icvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    genericRegionCoupledFluxFvPatchField<vector>(icvf, iF)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> regionCoupledTraction::flux() const
{
    // just a quick check that the used interface is a fsiInterface
    fsiInterface();

    // Lookup neighbouring patch field
    const volSymmTensorField& nbrSigmaField =
        nbrMesh().lookupObject<volSymmTensorField>
        (
            "sigma"
        );

    tmp<vectorField> tnbrTraction =
        nbrSigmaField.boundaryField()[neighbourPatchID()]
      & nbrPatch().nf();

    const vectorField& nbrTraction = tnbrTraction();

    // Calculate interpolated patch field
    vectorField tractionNbrToOwn = interpolateFromNbrField<vector>(nbrTraction);

    // Flip traction sign after transferring from fluid to solid
    tractionNbrToOwn = -tractionNbrToOwn;

    // Patch implicit stiffness field
    const fvPatchField<scalar>& impK =
        refPatch().lookupPatchField<volScalarField, scalar>("impK");

    // Patch reciprocal implicit stiffness field
    const fvPatchField<scalar>& rImpK =
        refPatch().lookupPatchField<volScalarField, scalar>("rImpK");

    // Patch gradient
    const fvPatchField<tensor>& gradD =
        refPatch().lookupPatchField<volTensorField, tensor>("gradD");

    // Patch Cauchy stress
    const fvPatchField<symmTensor>& sigma =
        refPatch().lookupPatchField<volSymmTensorField, symmTensor>("sigma");

    // Patch total deformation gradient inverse
    const fvPatchField<tensor>& Finv =
        refPatch().lookupPatchField<volTensorField, tensor>("Finv");

    //- Patch normals in initial/undeformed configuration
    const vectorField n(refPatch().nf());

    // Patch unit normals (deformed configuration)
    vectorField nCurrent(Finv.T() & n);
    nCurrent /= mag(nCurrent);

    // Tell the fsiInterface the total force
    fsiInterface().setTotalForce(gSum(tractionNbrToOwn*refPatch().magSf()));

    return
    (
        (
            tractionNbrToOwn
          - (nCurrent & sigma)
          + impK*(n & gradD)
        )*rImpK
    );
}

const regionInterfaces::fsiInterface&
regionCoupledTraction::fsiInterface() const
{
    if(   rgInterface().type()
       != regionInterfaces::fsiInterface::typeName )
    {
        FatalErrorInFunction
            << this->typeName << " BC can only "
            << "be used in combination with a "
            << regionInterfaces::fsiInterface::typeName
            << endl
            << exit(FatalError);
    }

    return refCast<const regionInterfaces::fsiInterface>
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
    regionCoupledTraction
);

} // End namespace Foam

// ************************************************************************* //
