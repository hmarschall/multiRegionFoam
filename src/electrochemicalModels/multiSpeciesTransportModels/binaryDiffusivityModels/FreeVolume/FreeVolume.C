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

#include "FreeVolume.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::FreeVolume<ThermoType>::FreeVolume
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
    rho2_("rho2", dimensionSet(1,-3,0,0,0,0,0), readScalar(binaryDiffDict_.lookup("soluteDensity"))),
    Kb_("Kb", dimensionSet(1,2,-2,-1,0,0,0), 1.380648e-23),
    Gamma_(0.8),
    unitsCorr_("units", dimensionSet(0,-0.5,0,0,0,0,0), 1)    
{ 

    //scalar Vf1_ = free volume[name2];
    scalar Vf1 = readScalar(binaryDiffDict_.subDict("freevolume").lookup(this->name1_));
    
    //scalar d1_ = diameter[name1];
    scalar d1 = readScalar(binaryDiffDict_.subDict("diameter").lookup(this->name1_));
    
   //scalar d2_ = diameter[name1];
    scalar d2 = readScalar(binaryDiffDict_.subDict("diameter").lookup(this->name2_));
    
     d2_= d2;
     A_= (Gamma_*d2)/(Vf1*d1);
     B_= d1/6;



}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::FreeVolume<ThermoType>::D() const
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

        d= B_*sqrt((3*Kb_*T)/(rho2_*pow(d2_,(3.0)))*exp(A_))*unitsCorr_;

    return tD;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
