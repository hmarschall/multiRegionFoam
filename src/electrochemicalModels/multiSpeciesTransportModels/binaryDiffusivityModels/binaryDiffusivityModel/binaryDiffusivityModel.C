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

#include "binaryDiffusivityModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// namespace Foam
// {
//     defineTypeNameAndDebug(binaryDiffusivityModel, 0);
//     defineRunTimeSelectionTable(binaryDiffusivityModel, dictionary);
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::binaryDiffusivityModel<ThermoType>::binaryDiffusivityModel
(
    const ThermoType& thermo,
    const word name1,
    const word name2,
    const IOdictionary& dict,
    const phaseModel& phase
)
:
    dict_(dict),
    name1_(name1),
    name2_(name2),
    // name1_(name1),
    // name2_(name2),
    thermo_(thermo),
    phase_(phase),
    phaseSys_(phase.fluid()),
    pZones_(phaseSys_.porosityModels())
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::autoPtr
<
    Foam::binaryDiffusivityModel<ThermoType>
>
Foam::binaryDiffusivityModel<ThermoType>::New
(
    const ThermoType& thermo,
    const word name1,
    const word name2,
    const IOdictionary& dict,
    const phaseModel& phase
)
{
    const word modelType
    (
        word(dict.lookup("binaryDiffusivityModel"))
        // + "<"
        // + word(ThermoType::typeName)
        // + ">"
    );

    Info<< "Selecting binaryDiffusivityModel " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown binaryDiffusivityModel type "
            << modelType << nl << nl
            << "Available types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<binaryDiffusivityModel>
    (
        cstrIter()(thermo, name1, name2, dict, phase)
    );
}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// template<class ThermoType>
// Foam::tmp<Foam::scalarField> Foam::binaryDiffusivityModel<ThermoType>::D
// (
//     const scalarField& p,
//     const scalarField& T,
//     const label patchi
// ) const
// {
//     notImplemented
//     (
//         "basicComponentThermo::D"
//         "(const scalarField& p, const scalarField& T, const label patchi) const"
//     );
//     return tmp<scalarField>(NULL);
// }

// template<class ThermoType>
// Foam::tmp<Foam::volScalarField> Foam::binaryDiffusivityModel<ThermoType>::D() const
// {
//     notImplemented("basicComponentThermo::D() const");
//     return volScalarField::null();
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

