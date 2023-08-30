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

#include "TynCalus.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::TynCalus<ThermoType>::TynCalus
(
    const ThermoType& thermo,
    const word& name1,
    const word& name2,
    const IOdictionary& dict,
    const phaseModel& phase
)
:
    binaryDiffusivityModel<ThermoType>(thermo, name1, name2, dict, phase),
    binaryDiffDict_(dict.subDict(typeName + "Coeff")),
    mu2_("mu2", dimensionSet(1,-1,-1,0,0,0,0), readScalar(binaryDiffDict_.lookup("solutetViscosity"))),
    unitsCorr_("units", dimensionSet(1, 1, -2, -1, 0, 0, 0), 1)
{   
    //scalar v1bp_ = MolarVolume[name1];
    //scalar v2bp_ = MolarVolume[name2];
    scalar V1bp = readScalar(binaryDiffDict_.subDict("MolarVolume").lookup(this->name1_));
    scalar V2bp = readScalar(binaryDiffDict_.subDict("MolarVolume").lookup(this->name2_));
    
    V12_ = (pow( V2bp, (0.267)) / (pow(V1bp, (0.43))));
   
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::TynCalus<ThermoType>::D() const
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

        d = 8.93e-12*V12_*T/mu2_* unitsCorr_;
    

    return tD;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
