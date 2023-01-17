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
    Foam::IQNILSRelaxation

Description
    Abstract base class for relaxation models

Author
    Constantin Habes, TU Darmstadt.  All rights reserved.

SourceFiles
    IQNILSRelaxation.C

\*---------------------------------------------------------------------------*/

#ifndef IQNILSRelaxation_H
#define IQNILSRelaxation_H

#include "relaxationModel.H"
#include "RectangularMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class IQNILSRelaxation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class IQNILSRelaxation
:
    public relaxationModel<Type>
{

private:

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const IQNILSRelaxation<Type>&);


protected:

    // Protected data

        //- Number of reused coupling data/ matrix columns
        //  from previous time steps
        const scalar reuse_; 

        //- Residual delta matrix
        DynamicList<Field<Type> > V_;

        //- BC value delta matrix
        DynamicList<Field<Type> > W_;

        //- List of time indices
        DynamicList<scalar> T_;

        //- Intermediate reference field for reaxed field
        Field<Type> fldRef_;

        //-Intermediate reference residual
        Field<Type> resRef_;



public:

    // Constructors

        //- Construct from dictionary
        IQNILSRelaxation
        (
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from runTime
        IQNILSRelaxation
        (
            const Time& runTime
        );

        //- Construct as copy
        IQNILSRelaxation
        (
            const IQNILSRelaxation<Type>&
        );


    // Destructor

        ~IQNILSRelaxation();


    // Member Functions

        // Edit
            //- Relax field 
            virtual void relax(Field<Type> &field);

            //- Update matrices V and W
            void updateVW(Field<Type> &field);

            //- Relax the field using IQN-ILS
            void relaxIQNILS(Field<Type> &field);



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "IQNILSRelaxation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
