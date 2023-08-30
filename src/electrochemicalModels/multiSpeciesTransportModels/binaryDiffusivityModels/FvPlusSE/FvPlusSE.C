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

#include "FvPlusSE.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::FvPlusSE<ThermoType>::FvPlusSE
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
    Kb_(1.380648e-23),
    mu_(readScalar(binaryDiffDict_.lookup("solventViscosity"))),
    rho1_(readScalar(binaryDiffDict_.lookup("solventDensity"))),
    rho2_(readScalar(binaryDiffDict_.lookup("soluteDensity"))),
    PI_(3.14159),
    Dself_Fv(readScalar(binaryDiffDict_.lookup("Dself_Fv"))),
    Dself_SE(readScalar(binaryDiffDict_.lookup("Dself_SE"))),
    unitsCorr_("units", dimensionSet(0, 2, -1, -1, 0, 0, 0), 1)
           
             
{  
    //scalar Vf1_ = free volume[name1];
    const scalar& Vf1 =readScalar(binaryDiffDict_.subDict("freevolume").lookup(this->name1_));    
   
    //scalar d1_ = diameter[name1];
    scalar d1 = readScalar(binaryDiffDict_.subDict("diameter").lookup(this->name1_));
    
   //scalar d2_ = diameter[name1];
    scalar d2 = readScalar(binaryDiffDict_.subDict("diameter").lookup(this->name2_));
   
  
    Gamma_ = 0.8;
    b_ = Dself_Fv/Dself_SE;    
    A_ = Kb_/3*PI_*mu_*d2 ;
    B_ = b_*(sqrt((d1*rho1_)/(d2*rho2_)));
    C_ = (-Gamma_/Vf1)*((d2/d1)-1);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::tmp<Foam::volScalarField> 
Foam::FvPlusSE<ThermoType>::D() const
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
    
	d = A_*T*(1+(B_*exp(C_)))* unitsCorr_;    

    return tD;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
