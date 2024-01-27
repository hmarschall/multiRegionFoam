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

#include "Ostream.H"
#include "error.H"
#include "regionCoupledPressureValue.H"
#include "addToRunTimeSelectionTable.H"
#include "regionInterfaceType.H"
#include "capillaryInterface.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupledPressureValue::
regionCoupledPressureValue
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    genericRegionCoupledJumpFvPatchField<scalar>(p, iF)
{}


Foam::regionCoupledPressureValue::
regionCoupledPressureValue
(
    const regionCoupledPressureValue& icpv,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    genericRegionCoupledJumpFvPatchField<scalar>(icpv, p, iF, mapper)
{}


Foam::regionCoupledPressureValue::
regionCoupledPressureValue
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    genericRegionCoupledJumpFvPatchField<scalar>(p, iF, dict)
{}


Foam::regionCoupledPressureValue::
regionCoupledPressureValue
(
    const regionCoupledPressureValue& icpv,
    const DimensionedField<scalar, volMesh>& iF
)
:
    genericRegionCoupledJumpFvPatchField<scalar>(icpv, iF)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<scalarField> regionCoupledPressureValue::valueJump() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // surface tension
    const areaScalarField& sigma = capInterface().sigma();

    // interfacial curvature
    const areaScalarField& K = capInterface().aMesh().faceCurvatures();


    // surface velocity terms
    const areaVectorField& Us = capInterface().Us();

    areaScalarField divUs(fac::div(Us));
    divUs.correctBoundaryConditions();

    // gravity term
    vector pRefPoint(mesh.solutionDict().subDict("PISO").lookup("pRefPoint"));

    dimensionedVector g (capInterface().gravitationalProperties().lookup("g"));

    dimensionedScalar muFluidNbr
    (
        nbrMesh().lookupObject<IOdictionary>("transportProperties")
        .lookup("mu")
    );

    dimensionedScalar muFluid
    (
        mesh.lookupObject<IOdictionary>("transportProperties")
        .lookup("mu")
    );

    dimensionedScalar rhoFluidNbr
    (
        nbrMesh().lookupObject<IOdictionary>("transportProperties")
        .lookup("rho")
    );

    dimensionedScalar rhoFluid
    (
        mesh.lookupObject<IOdictionary>("transportProperties")
        .lookup("rho")
    );

    return
    (
        2.0*(muFluidNbr.value() - muFluid.value())*divUs.internalField()
      - sigma.internalField()*K.internalField()
      + (rhoFluidNbr.value() - rhoFluid.value())
        *(
            (
                mesh.Cf().boundaryField()[this->patch().index()]
              - pRefPoint
            ) & g.value()
        )
   );
}

const regionInterfaces::capillaryInterface&
regionCoupledPressureValue::capInterface() const
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
    regionCoupledPressureValue
);

} // End namespace Foam

// ************************************************************************* //
