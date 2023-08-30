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

#include "multiSpeciesTransportModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::multiSpeciesTransportModel>
Foam::multiSpeciesTransportModel::New
(
    const phaseModel& phase
)
{
    word multiSpeciesTransportModelType;
    // word multiSpeciesTransportModelType = "MaxwellStefan<" + phase.thermo().type() + ">";

        // Enclose the creation of the dictionary to ensure it is deleted
    // before the multiSpeciesTransportModel is created otherwise the
    // dictionary is entered in the database twice
    {
        IOdictionary dict
        (
            IOobject
            (
                IOobject::groupName("thermophysicalTransport", phase.name()),
                phase.mesh().time().constant(),
                phase.mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        multiSpeciesTransportModelType
        =
        (
            word(dict.lookup("type"))
            // + "<"
            // + phase.thermo().type()
            // + ">"
        );
    }

    Info<< "Selecting multiSpeciesTransportModel: "
    << multiSpeciesTransportModelType << endl;
  
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(multiSpeciesTransportModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "MultiSpeciesDiffusivityModel::New(const volVectorField&, "
            "const surfaceScalarField&)"
        )   << "Unknown multiSpeciesDiffusivityModel type "
            << multiSpeciesTransportModelType << endl << endl
            << "Valid  multiSpeciesDiffusivityModel are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
  }

  return autoPtr<multiSpeciesTransportModel>
      (cstrIter()(phase));
}


// ************************************************************************* //
