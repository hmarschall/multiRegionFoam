/*---------------------------------------------------------------------------*\u87
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

#include "Fuller.H"
// #include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::Fuller<ThermoType>::Fuller
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
    unitsCorr_("units", dimensionSet(1, 1, -3, -1.75, 0, 0, 0), 1)
{   

    label speciesIndex1 = this->thermo_.composition().species()[this->name1_];
    label speciesIndex2 = this->thermo_.composition().species()[this->name2_];

    const scalar& W1 = this->thermo_.composition().Wi(speciesIndex1);
    const scalar& W2 = this->thermo_.composition().Wi(speciesIndex2);

    // TODO: Those values also have to be automaticially provided by the solver itself
    // const scalar& V1 = readScalar(binaryDiffDict_.subDict("diffusionVolume").lookup(this->name1_));
    // const scalar& V2 = readScalar(binaryDiffDict_.subDict("diffusionVolume").lookup(this->name2_));
    const scalar& V1 = Foam::diffusivityModels::fsgDiffusionVolumes(this->name1_);
    const scalar& V2 = Foam::diffusivityModels::fsgDiffusionVolumes(this->name2_);

    W12 = (W1 * W2) / (W1 + W2);
    V12 = sqr( pow( V1, (1.0/3.0) ) + pow( V2, (1.0/3.0)) );
    pAtm =  101325;
}

// ====================== Fuller ========================== //
//             1e-3 * T^{1.75} * sqrt(1/mA + 1/mB)
//  D = 1e-4 * -----------------------------------
//                p * [ vA^(1/3) + vB^{1/3} ]^2
//  where
//      D = diffusivity ......... [m^2/s]
//      T = temperature ......... [K]
//      p = total pressure ...... [atm]
//      m = molecular weight .... [kg/kmol]
//      v = diffusion volume .... [cm^3]     NOTE: cm
//      A,B = species index
//
//  Fuller, Schettler, and Giddings,
//  A new method for prediction of binary gas-phase diffusion coefficients,
//  Industrial and Engineering Chemistry, v58, n5, May, 1966, pp 19-27.



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::Fuller<ThermoType>::D() const
{
    const volScalarField& T = this->thermo_.T();
    const volScalarField& p = this->thermo_.p();

    const fvMesh& mesh = T.mesh();

    // Info << "Print this->name1_ " << this->name1_ << endl;
    // Info << "Print this->name2_ " << this->name2_ << endl;

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

    volScalarField& d = tD();

    volScalarField pTot = p/pAtm;

    // TODO: - Get rid of these forAll loop over the cells
    // Due to the lacking of unit conservative:
    // Keep units of the fields (locally varying variables)
    // All constants are dimensionless
    // Compensate the missing dimensions to the wished m^2/s via a dimensionedScalar

    d = (1.011e-7 * pow(T,(1.75)) / sqrt(W12) / pTot / V12) * unitsCorr_;

    Info << "Print D_ " + this->name1_ + "_" + this->name2_ << endl;
    Info << "Print average(D_) " << average(d) << endl;
    // Info << "Print d " << d << endl;

    return tD;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
