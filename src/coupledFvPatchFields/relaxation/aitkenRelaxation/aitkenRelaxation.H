/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  Field<Type> ield         | foam-extend: Open Source CFD
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

Class
    Foam::aitkenRelaxation

Description
    Abstract base class for relaxation models

Author
    Constantin Habes, TU Darmstadt.  All rights reserved.

SourceFiles
    aitkenRelaxation.C

\*---------------------------------------------------------------------------*/

#ifndef fixedRelaxation_H
#define fixedRelaxation_H

#include "relaxationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class aitkenRelaxation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class aitkenRelaxation
:
    public relaxationModel<Type>
{

private:

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const aitkenRelaxation<Type>&);


protected:

    // Protected data

        //- Aitken relaxation factor
        scalar aitkenRelax_;

public:

    // Constructors

        //- Construct from dictionary
        aitkenRelaxation
        (
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from runTime
        aitkenRelaxation
        (
            const Time& runTime
        );

        //- Construct as copy
        aitkenRelaxation
        (
            const aitkenRelaxation<Type>&
        );


    // Destructor

        ~aitkenRelaxation();


    // Member Functions

        // Edit
            //- Relax field 
            virtual void relax(Field<Type> &field);

            //- Update Aitken relaxation factor
            void updateAitkenFactor(Field<Type> &field);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "aitkenRelaxation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
