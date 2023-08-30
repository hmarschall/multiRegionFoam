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

#include "StokeEinstein.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::StokeEinstein<ThermoType>::StokeEinstein
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
    Kb_("Kb", dimensionSet(1,2,-2,-1,0,0,0), 1.380648e-23),
    mu_("mu", dimensionSet(1,-1,-1,0,0,0,0), readScalar(binaryDiffDict_.lookup("solventViscosity"))),
    d_("d", dimensionSet(1,2,-2,-1,0,0,0), readScalar(binaryDiffDict_.lookup("sphericalDiameter"))),
    PI_(3.14159),
    unitsCorr_("units", dimensionSet(-3, -1, 4, 1, 0, 0, 0), 1)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::StokeEinstein<ThermoType>::D() const
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

    // Unit conservative [m^2/s]
    d = Kb_*T/3*PI_*mu_*d_*unitsCorr_;

    return tD;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
