/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "phaseModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseModel, 0);
    defineRunTimeSelectionTable(phaseModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseModel::phaseModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    dict_(dict),
    volScalarField
    (
        IOobject
        (
            "alpha",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh()
    ),
    porosityModels_
    (
       new porosityModelList
       (
            mesh
       ) 
    ),
    name_(mesh.name())
{

}


Foam::autoPtr<Foam::phaseModel> Foam::phaseModel::clone() const
{
    NotImplemented;
    return autoPtr<phaseModel>(nullptr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseModel::~phaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvMesh& Foam::phaseModel::mesh() const
{
    return mesh_;
}

const Foam::dictionary& Foam::phaseModel::dict() const
{
    return dict_;
}

const Foam::word& Foam::phaseModel::name() const
{
    return name_;
}

void Foam::phaseModel::correct()
{
}


void Foam::phaseModel::correctKinematics()
{}


void Foam::phaseModel::correctThermo()
{}


void Foam::phaseModel::correctTurbulence()
{}


void Foam::phaseModel::correctEnergyTransport()
{}


// ************************************************************************* //
