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
    Foam::fixedRelaxation

Description
    Abstract base class for relaxation models

Author
    Constantin Habes, TU Darmstadt.  All rights reserved.

SourceFiles
    fixedRelaxation.C

\*---------------------------------------------------------------------------*/

#ifndef fixedRelaxation_H
#define fixedRelaxation_H

#include "relaxationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fixedRelaxation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixedRelaxation
:
    public relaxationModel<Type>
{
protected:

    // Protected data

        //- Fixed relaxation factor
        scalar relax_;

public:

    //- Runtime type information
    TypeName("fixed");

    // Constructors

        //- Construct from dictionary
        fixedRelaxation
        (
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct as copy
        fixedRelaxation
        (
            const fixedRelaxation<Type>&
        );

        virtual autoPtr<relaxationModel<Type> > clone() const
        {
            return autoPtr<relaxationModel<Type> >
            (
                new fixedRelaxation<Type>(*this)
            );
        }

    // Destructor

        ~fixedRelaxation();


    // Member Functions

        // Edit
            virtual void initialize(const Field<Type> &field);
            
            //- Relax field 
            virtual void relax(Field<Type> &field);

        //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fixedRelaxation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
