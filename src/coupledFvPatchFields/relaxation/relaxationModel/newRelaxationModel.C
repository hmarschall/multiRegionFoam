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
#include "relaxationModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::relaxationModel<Type> >
Foam::relaxationModel<Type>::New
(
    const Time& runTime,
    const dictionary& dict
)
{
    if (dict.empty())
    {
            WarningIn
            (
                "relaxationModel<Type>::New\n"
            )
                << "relaxationModel is constructed without dictionary.\n"
                << "\tSelecting a fixed relaxation model with relaxation factor of 1 "
                << "to be safe."
                << endl;
    } 
    else if (!dict.found("relaxType")) 
    {
            WarningIn
            (
                "relaxationModel<Type>::New\n"
            )
                << "No relaxType entry in " << dict.dictName() << "\n"
                << "\tAssuming a fixed relaxation model with relaxation factor of 1 "
                << endl;
    }

    word relaxationModelTypeName
        (
            dict.lookupOrDefault<word>("relaxType", "fixed")
        );

    if (debug)
    {
        Info<< "relaxationModel<Type>::New(const Time& runTime, "
               "const dictionary& dict) : relaxationModelTypeName = "  
            << relaxationModelTypeName
            << endl;
    }

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(relaxationModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "relaxationModel::New(const Time& runTime, const dictionary& dict) "
        )   << "Unknown relaxationModel type "
            << relaxationModelTypeName << nl << nl
            << "Valid relaxationModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<relaxationModel<Type> >(cstrIter()(runTime, dict));
}

// ************************************************************************* //