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
    Foam::relaxationModel

Description
    Abstract base class for relaxation models

Author
    Constantin Habes, TU Darmstadt.  All rights reserved.

SourceFiles
    relaxationModel.C

\*---------------------------------------------------------------------------*/

#ifndef relaxationModel_H
#define relaxationModel_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "foamTime.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class relaxationModel Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class relaxationModel
{
protected:

    // Protected data
    
        const Time& runTime_;

        //- Curent time
        scalar curTime_;

        //- Counter for applied correction steps 
        label corr_;

        //- Previous field
        Field<Type> prevFld_;

        //- Residual field
        Field<Type> resFld_;

        //- Previous residual field
        Field<Type> prevResFld_;

        //- Initial relaxation factor
        const scalar initRelax_;

public:

    //- Runtime type information
    TypeName("relaxationModel");


    //- Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            relaxationModel,
            dictionary,
            (
                const Time& runTime,
                const dictionary& dict
            ),
            (runTime, dict)
        );


    // Constructors

        //- Construct from dictionary
        relaxationModel
        (
            const Time& runTime,
            const dictionary& dict = dictionary()
        );

        //- Construct as copy
        relaxationModel
        (
            const relaxationModel<Type>&
        );

        //- Construct and return a clone
        virtual autoPtr<relaxationModel<Type> > clone() const = 0;

    // Selectors

        //- Return a pointer to the selected relaxation model
        //  created from a dictionary
        static autoPtr<relaxationModel> New
        (
            const Time& runTime,
            const dictionary& dict = dictionary()
        );


    // Destructor

        virtual ~relaxationModel();


    // Member Functions

        // Access

            //- Return const current residual
            const Field<Type>& residual()
            {
                return resFld_;
            };

            //- Return const initial relaxation factor
            scalar initRelax()
            {
                return initRelax_;
            };

        // Edit
            virtual void initialize(const Field<Type> &field);

            //- Relax field 
            virtual void relax(Field<Type> &field) = 0;

            //- Compute residual
            void updateResiual(const Field<Type> &field);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeRelaxationModel(Type)                                             \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(relaxationModel<Type>, 0);            \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        relaxationModel<Type>,                                                \
        dictionary                                                            \
    );


#define makeRelaxationModelType(RS, Type)                                     \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(RS<Type>, 0);                         \
                                                                              \
    relaxationModel<Type>::adddictionaryConstructorToTable<RS<Type> >       \
        add##RS##Type##ConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "relaxationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
