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
    An acceleration model that applies a dynamic relaxation factor based on
    the Aitken procedure

Author
    Constantin Habes, TU Darmstadt.  All rights reserved.

SourceFiles
    aitkenRelaxation.C

\*---------------------------------------------------------------------------*/

#ifndef aitkenRelaxation_H
#define aitkenRelaxation_H

#include "accelerationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class aitkenRelaxation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class aitkenRelaxation
:
    public accelerationModel<Type>
{
protected:

    // Protected data

        //- Aitken relaxation factor
        scalar aitkenRelax_;

public:

    //- Runtime type information
    TypeName("aitken");

    // Constructors

        //- Construct from dictionary
        aitkenRelaxation
        (
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct as copy
        aitkenRelaxation
        (
            const aitkenRelaxation<Type>&
        );

        //- Construct and return a clone
        virtual autoPtr<accelerationModel<Type> > clone() const
        {
            return autoPtr<accelerationModel<Type> >
            (
                new aitkenRelaxation<Type>(*this)
            );
        }


    // Destructor

        ~aitkenRelaxation();


    // Member Functions

        // Edit
            virtual void initialize(const Field<Type> &field);

            //- Relax field
            virtual void relax(Field<Type> &field);

            //- Update Aitken relaxation factor
            void updateAitkenFactor(Field<Type> &field);

        //- Write
            virtual void write(Ostream&) const;

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
