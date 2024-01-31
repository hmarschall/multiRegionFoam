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

#include "error.H"
#include "accelerationModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::accelerationModel<Type> >
Foam::accelerationModel<Type>::New
(
    const Time& runTime,
    const dictionary& dict
)
{
    if (dict.empty())
    {
            WarningIn
            (
                "accelerationModel<Type>::New\n"
            )
                << "accelerationModel is constructed without dictionary.\n"
                << "\tSelecting a fixed relaxation model with a "
                << "relaxation factor of 1 "
                << "to be safe."
                << endl;
    }
    else if (!dict.found("accType"))
    {
            WarningIn
            (
                "accelerationModel<Type>::New\n"
            )
                << "No accType entry in " << dict.dictName() << "\n"
                << "\tAssuming that no acceleration model should be applied"
                << endl;
    }

    word accelerationModelTypeName
        (
            dict.lookupOrDefault<word>("accType", "fixed")
        );

    if (debug)
    {
        Info<< "accelerationModel<Type>::New(const Time& runTime, "
               "const dictionary& dict) : accelerationModelTypeName = "
            << accelerationModelTypeName
            << endl;
    }

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(accelerationModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "accelerationModel::New(const Time& runTime, const dictionary& dict)"
        )   << "Unknown accelerationModel type "
            << accelerationModelTypeName << nl << nl
            << "Valid accelerationModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<accelerationModel<Type> >(cstrIter()(runTime, dict));
}

// ************************************************************************* //