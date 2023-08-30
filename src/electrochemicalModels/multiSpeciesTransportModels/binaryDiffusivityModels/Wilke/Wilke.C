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

#include "Wilke.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::Wilke<ThermoType>::Wilke
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
    T0_("T0", dimTemperature, 1),
    unitsCorr_("units", dimensionSet(1, 1, -3, -1.5, 0, 0, 0), 1)
{

    A = 1.06036;    B = 0.15610;    C = 0.19300;    DD = 0.47635;
    E = 1.03587;    F = 1.52996;    G = 1.76474;    H = 3.89411;   

    const scalar& epsLJ1 = readScalar(binaryDiffDict_.subDict("epsilonLJ").lookup(this->name1_));
    const scalar& epsLJ2 = readScalar(binaryDiffDict_.subDict("epsilonLJ").lookup(this->name2_));

    const scalar& sigma1 =
        readScalar(binaryDiffDict_.subDict("collisionalDiametre").lookup(this->name1_));
    const scalar& sigma2 =
        readScalar(binaryDiffDict_.subDict("collisionalDiametre").lookup(this->name2_));

    label speciesIndex1 = this->thermo_.composition().species()[this->name1_];
    label speciesIndex2 = this->thermo_.composition().species()[this->name2_];

    const scalar& W1 = this->thermo_.composition().W(speciesIndex1);
    const scalar& W2 = this->thermo_.composition().W(speciesIndex2);

    sigma_ij = (sigma1 + sigma2) / 2;
    W12 = (W1 * W2) / (W1 + W2);
    phi =  101325;
    E_ij = sqrt(epsLJ1*epsLJ2);

}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::Wilke<ThermoType>::D() const
{
    const volScalarField& T = this->thermo_.T();
    const volScalarField& p = this->thermo_.p();

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

     tmp<volScalarField> T_ND
    (
        new volScalarField
        (
            IOobject
            (
                "TN",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimless
        )
    );

    // dimensionless collision integral
    tmp<volScalarField> omegaD
    (
        new volScalarField
        (
            IOobject
            (
                "omega",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimless
        )
    );

    volScalarField& d = tD.ref();
    volScalarField& T_N = T_ND.ref();
    volScalarField& omega = omegaD.ref();

    T_N = (T/T0_) / E_ij;

    omega = A/pow(T_N,B) + C/exp(DD*T_N) + E/exp(F*T_N) + G/exp(H*T_N);
    
    d = (1.858e-7 * phi * sqrt(pow(T,3) / W12)
            / p / sqr(sigma_ij) / omega) * unitsCorr_;


    return tD;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
