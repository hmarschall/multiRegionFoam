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

#include "regionCoupledHeatFlux.H"
#include "regionCoupledScalarJump.H"
#include "addToRunTimeSelectionTable.H"
#include "heatTransferInterface.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupledHeatFlux::
regionCoupledHeatFlux
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    genericRegionCoupledFluxFvPatchField<scalar>(p, iF)
{}


Foam::regionCoupledHeatFlux::
regionCoupledHeatFlux
(
    const regionCoupledHeatFlux& icpf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    genericRegionCoupledFluxFvPatchField<scalar>(icpf, p, iF, mapper)
{}


Foam::regionCoupledHeatFlux::
regionCoupledHeatFlux
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    genericRegionCoupledFluxFvPatchField<scalar>(p, iF, dict)
{}


Foam::regionCoupledHeatFlux::
regionCoupledHeatFlux
(
    const regionCoupledHeatFlux& icpf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    genericRegionCoupledFluxFvPatchField<scalar>(icpf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> regionCoupledHeatFlux::fluxJump() const
{
    return (*this * 0);
}

const regionInterfaces::heatTransferInterface&
regionCoupledHeatFlux::heatTransInterface() const
{
    if(   rgInterface().type()
       != regionInterfaces::heatTransferInterface::typeName )
    {
        FatalErrorInFunction
            << this->typeName << " BC can only "
            << "be used in combination with a "
            << regionInterfaces::heatTransferInterface::typeName
            << endl
            << exit(FatalError);
    }

    return refCast<const regionInterfaces::heatTransferInterface>
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
    regionCoupledHeatFlux
);

} // End namespace Foam

// ************************************************************************* //
