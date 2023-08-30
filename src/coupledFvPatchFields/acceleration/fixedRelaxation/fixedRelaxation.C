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

\*---------------------------------------------------------------------------*/

#include "fixedRelaxation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class Type>
Foam::fixedRelaxation<Type>::fixedRelaxation
(
    const Time& runTime,
    const dictionary& dict
)
:
    accelerationModel<Type>(runTime, dict),
    relax_(this->initRelax_)
{
    if (relax_ != 1.0)
    {
        Info<< "Selecting an fixedRelaxation model for " << dict.dictName()
            << " with fixed relaxation factor " << this->relax_
            << endl;
    }
}

template<class Type>
Foam::fixedRelaxation<Type>::fixedRelaxation
(
    const fixedRelaxation<Type>& fR
)
:
    accelerationModel<Type>(fR),
    relax_(fR.relax_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::fixedRelaxation<Type>::~fixedRelaxation()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::fixedRelaxation<Type>::initialize(const Field<Type> &curFld)
{
    Foam::accelerationModel<Type>::initialize(curFld);
}

template<class Type>
void Foam::fixedRelaxation<Type>::relax(Field<Type> &curFld)
{
    //- Check if solver time is time saved by relaxation model

    if (this->runTime_.value() != this->curTime_)
    {
        //- Reset counter for corrector steps
        this->corr_ = 1;
        //- Set curent time to solver time
        this->curTime_ = this->runTime_.value();
    }

    //- Update residuals
    this->updateResiual(curFld);

    if (relax_ != 1.0)
    {
        Info<< nl
            << "Relaxing field with fixed relaxation factor: "
            << relax_
            << endl;
    }

    //- Relax field
    curFld = this->prevFld_ + this->initRelax_ * this->resFld_;

    //- Store relaxed field a new field
    this->prevFld_ = curFld;

    //- Increment corrector step counter
    this->corr_++;
}

template<class Type>
void Foam::fixedRelaxation<Type>::write(Ostream& os) const
{
    accelerationModel<Type>::write(os);
}

// ************************************************************************* //
