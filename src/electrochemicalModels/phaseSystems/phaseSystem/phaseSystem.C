/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "phaseSystem.H"
// #include "fvcDiv.H"
// #include "fvcGrad.H"
// #include "fvcSnGrad.H"
#include "pimpleControl.H"
// #include "pressureReference.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseSystem, 0);
    defineRunTimeSelectionTable(phaseSystem, dictionary);
}

const Foam::word Foam::phaseSystem::propertiesName("phaseProperties");


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseSystem::phaseSystem
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            propertiesName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    mesh_(mesh),
    phaseModels_(0),
    // phaseModels_
    // (
    //     lookup("phases"),
    //     phaseModel::iNew(*this)
    // ),
    porosityModels_
    (
       new porosityModelList
       (
            mesh
       ) 
    ),
    p_rgh_
    (
        IOobject
        (
            "p_rgh",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    pRefCell_(0),
    pRefValue_(0.0),
    pimple_(const_cast<fvMesh&>(mesh))
{
    reset(*this);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseSystem::~phaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseSystem::reset(const IOdictionary& dict)
{
    wordList phaseNames;

    label j = 0;

    HashTable<wordList> phaseList = dict.lookup("phaseNames");

    forAllConstIter(HashTable<wordList>, phaseList, iter)
    {
        const wordList& phases = iter();

        forAll(phases, phaseI)
        {

            if (findIndex(phaseNames, phases[phaseI]))
            {
                phaseNames.setSize(phaseNames.size()+1);
                phaseNames[j] = phases[phaseI];
            }

            j++;
        }
    }

    phaseModels_.setSize(phaseNames.size());

    label i = 0;

    forAllConstIter(HashTable<wordList>, phaseList, iter)
    {
        // const word& modelType = iter.key();
        const wordList& phases = iter();

        if (phases.size())
        {
            forAll(phases, phaseI)
            {
                Info << "Creating " << phases[phaseI] << endl;

                phaseModels_.set
                (
                    i++,
                    phaseModel::New
                    (
                        *this,
                        phases[phaseI],
                        i-1
                    )
                );
            }
        }
    }
}

void Foam::phaseSystem::correct()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correct();
    }
}

void Foam::phaseSystem::correctThermo()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctThermo();
    }
}

void Foam::phaseSystem::correctEnergyTransport()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctEnergyTransport();
    }
}


// ************************************************************************* //
