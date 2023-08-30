/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "multiComponentMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::constructSpeciesData
(
    const dictionary& thermoDict
)
{
    forAll(species_, i)
    {
        speciesData_.set
        (
            i,
            new ThermoType(thermoDict.subDict(species_[i]))
        );
    }

    // Info << "Debug constructSpeciesData in multiComponentMixture " << endl;

    return speciesData_[0];
}


template<class ThermoType>
void Foam::multiComponentMixture<ThermoType>::correctMassFractions()
{
    // Multiplication by 1.0 changes Yt patches to "calculated"
    volScalarField Yt("Yt", 1.0*Y_[0]);

    for (label n=1; n<Y_.size(); n++)
    {
        Yt += Y_[n];
    }

    if (mag(max(Yt).value()) < ROOTVSMALL)
    {
        FatalErrorInFunction
            << "Sum of mass fractions is zero for species " << this->species()
            << exit(FatalError);
    }

    forAll(Y_, n)
    {
        Y_[n] /= Yt;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::multiComponentMixture<ThermoType>::multiComponentMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const HashPtrTable<ThermoType>& specieThermoData,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecieMixture(thermoDict, specieNames, mesh, phaseName),
    speciesData_(species_.size()),
    mixture_("mixture", *specieThermoData[specieNames[0]]),
    mixtureVol_("volMixture", *specieThermoData[specieNames[0]])
{
    // Set the specie template type, e.g. sutherlandTransport<specieThermo<janafThermo<perfectGas> > >
    forAll(species_, i)
    {
        speciesData_.set
        (
            i,
            new ThermoType(*specieThermoData[species_[i]])
        );
    }

    correctMassFractions();
}


template<class ThermoType>
Foam::multiComponentMixture<ThermoType>::multiComponentMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecieMixture(thermoDict, thermoDict.lookup("species"), mesh, phaseName),
    speciesData_(species_.size()),
    mixture_("mixture", constructSpeciesData(thermoDict)),
    mixtureVol_("volMixture", speciesData_[0])
{
    // // Info << "Print constructor multiComponentMixture 1" << endl;
    correctMassFractions();
    // // Info << "Print constructor multiComponentMixture 2" << endl;
}


// // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::cellMixture
(
    const label celli
) const
{
    // mixture_ = Y_[0][celli]/speciesData_[0].W()*speciesData_[0];

    // // What is the purpose of this here
    // // ->> uses somehow the overloaded += operator for Thermotype mixture_, but dunno for what
    // // Info << "Print multiComp cellMixture " << endl;
    // // Info << "Print Y_.size() " << Y_.size() << endl;

    // for (label n=1; n<Y_.size(); n++)
    // {
    //     Info << "Print multiComp cellMixture davor " << endl;
    //     // Info << "Print multiComp speciesData_[n].W() " << speciesData_[n].W() << endl;
    //     // Info << "Print multiComp speciesData_[n] " << speciesData_[n] << endl;
    //     // Info << "Print multiComp Y_[n][celli]/speciesData_[n].W()*speciesData_[n] " << Y_[n][celli]/speciesData_[n].W()*speciesData_[n] << endl;
    //     mixture_ += Y_[n][celli]/speciesData_[n].W()*speciesData_[n];
    //     Info << "Print multiComp cellMixture danach " << endl;
    //     // Info << "Print Y_[n].name() " << Y_[n].name() << endl;
    //     // Info << "Print multiComp speciesData_[n].W() " << speciesData_[n].W() <<  endl;
    // }

    // // Info << "Print multiComp cellMixture 1" << endl;

    // return mixture_;

    mixture_ = Y_[0][celli]*speciesData_[0];

    // Info << "Print multiComp Y_[0][celli] " << Y_[0][celli] <<  endl;
    // Info << "Print multiComp cellMixture 2" << endl;

    for (label n=1; n<Y_.size(); n++)
    {
        // Info << "Print Y_[n].name() " << Y_[n].name() << endl;
        // Info << "Print Y_[n][celli] " << Y_[n][celli] << endl;
        // Info << "Print speciesData_[n] " << speciesData_[n] << endl;
        mixture_ += Y_[n][celli]*speciesData_[n];
        // Info << "Print Y_[n].name() " << Y_[n].name() << endl;
    }

    // Info << "Print multiComp cellMixture 3" << endl;

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::patchFaceMixture
(
    const label patchi,
    const label facei
) const
{
    // mixture_ =
    //     Y_[0].boundaryField()[patchi][facei]
    //    /speciesData_[0].W()*speciesData_[0];

    // for (label n=1; n<Y_.size(); n++)
    // {
    //     mixture_ +=
    //         Y_[n].boundaryField()[patchi][facei]
    //        /speciesData_[n].W()*speciesData_[n];
    // }

    // return mixture_;
    // Info << "Print Y_[0].boundaryField()[patchi][facei] " << Y_[0].boundaryField()[patchi][facei] << endl;
    mixture_ = Y_[0].boundaryField()[patchi][facei]*speciesData_[0];
    // Info << "Print multiComp patchFaceMixture 1" << endl;

    for (label n=1; n<Y_.size(); n++)
    {
        mixture_ += Y_[n].boundaryField()[patchi][facei]*speciesData_[n];
    }
    // Info << "Print multiComp patchFaceMixture 2" << endl;

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::cellVolMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    scalar rhoInv = 0.0;
    forAll(speciesData_, i)
    {
        rhoInv += Y_[i][celli]/speciesData_[i].rho(p, T);
    }

    mixtureVol_ =
        Y_[0][celli]/speciesData_[0].rho(p, T)/rhoInv*speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixtureVol_ +=
            Y_[n][celli]/speciesData_[n].rho(p, T)/rhoInv*speciesData_[n];
    }

    return mixtureVol_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::
patchFaceVolMixture
(
    const scalar p,
    const scalar T,
    const label patchi,
    const label facei
) const
{
    scalar rhoInv = 0.0;
    forAll(speciesData_, i)
    {
        rhoInv +=
            Y_[i].boundaryField()[patchi][facei]/speciesData_[i].rho(p, T);
    }

    mixtureVol_ =
        Y_[0].boundaryField()[patchi][facei]/speciesData_[0].rho(p, T)/rhoInv
      * speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixtureVol_ +=
            Y_[n].boundaryField()[patchi][facei]/speciesData_[n].rho(p,T)
          / rhoInv*speciesData_[n];
    }

    return mixtureVol_;
}



template<class ThermoType>
void Foam::multiComponentMixture<ThermoType>::read
(
    const dictionary& thermoDict
)
{
    forAll(species_, i)
    {
        // speciesData_[i] = ThermoType(thermoDict.lookup(species_[i]));
        speciesData_[i] = ThermoType(thermoDict.subDict(species_[i]));
    }
}


// ************************************************************************* //
