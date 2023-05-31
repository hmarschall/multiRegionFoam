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

#include "aitkenRelaxation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class Type>
Foam::aitkenRelaxation<Type>::aitkenRelaxation
(
    const Time& runTime,
    const dictionary& dict
)
:
    accelerationModel<Type>(runTime, dict),
    aitkenRelax_(this->initRelax_)
{
    Info<< "Selecting an aitkenRelaxation model for " << dict.dictName() 
        << " with initial relaxation factor " << this->initRelax_
        << endl;
}

template<class Type>
Foam::aitkenRelaxation<Type>::aitkenRelaxation
(
    const aitkenRelaxation<Type>& fR
)
:
    accelerationModel<Type>(fR),
    aitkenRelax_(fR.aitkenRelax_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::aitkenRelaxation<Type>::~aitkenRelaxation()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


template<class Type>
void Foam::aitkenRelaxation<Type>::initialize(const Field<Type> &curFld)
{
    Foam::accelerationModel<Type>::initialize(curFld);
}

template<class Type>
void Foam::aitkenRelaxation<Type>::relax(Field<Type> &curFld)
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

    //- Update Aitken relaxation factor
    updateAitkenFactor(curFld);

    Info<< nl
        << "Relaxing field with Aitken relaxation factor: "
        << aitkenRelax_ 
        << endl;

    //- Relax field
    curFld = this->prevFld_ + aitkenRelax_ * this->resFld_;

    //- Store relaxed field as new field
    this->prevFld_ = curFld;

    //- Increment corrector step counter
    this->corr_++;
}

template<class Type>
void Foam::aitkenRelaxation<Type>::write(Ostream& os) const
{
    accelerationModel<Type>::write(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "aitkenRelaxationTemplates.C"
#endif

// ************************************************************************* //
