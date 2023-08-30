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

#include "WilkeChang.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::WilkeChang<ThermoType>::WilkeChang
(
    const ThermoType& thermo,
    const word name1,
    const word name2,
    const IOdictionary& dict,
    const phaseModel& phase
)
:
    binaryDiffusivityModel<ThermoType>(thermo, name1, name2, dict, phase),
    binaryDiffDict_(dict.subDict(typeName + "Coeff")),
    v1_(readScalar(binaryDiffDict_.lookup("soluteMolarVolume"))),
    phi2_(readScalar(binaryDiffDict_.lookup("associationFactor"))),
    mu2_(readScalar(binaryDiffDict_.lookup("solventDynamicViscosity"))),
    W2_(readScalar(binaryDiffDict_.lookup("solventMolarWeigth"))),
    unitsCorr_("units", dimensionSet(0, 2, -1, -1, 0, 0, 0), 1)
{

    // label speciesIndex1 = this->thermo_.composition().species()[name1];
    // label speciesIndex2 = this->thermo_.composition().species()[name2];

    // const scalar& W1_ = this->thermo_.composition().Wi(speciesIndex1);
    // const scalar& W2_ = this->thermo_.composition().Wi(speciesIndex2);

}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::WilkeChang<ThermoType>::D() const
{
    const volScalarField& T = this->thermo_.T();

    const fvMesh& mesh = T.mesh();

    tmp<volScalarField> tD
    (
        new volScalarField
        (
            IOobject
            (
                "D_" + this->name1_ + "_" + this->name2_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -1, 0, 0)
        )
    );

    volScalarField& d = tD.ref();

    d =  (T*7.4e-12*sqrt(phi2_*W2_)/mu2_*pow(v1_,(0.6))) * unitsCorr_;

    return tD;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
