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

#include "IQNILSRelaxation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class Type>
Foam::IQNILSRelaxation<Type>::IQNILSRelaxation
(
    const Time& runTime,
    const dictionary& dict
)
:
    relaxationModel<Type>(runTime, dict),
    reuse_(dict.lookupOrDefault<label>("couplingReuse", 0)),
    V_(),
    W_(),
    T_(),
    fldRef_(),
    resRef_()
{
    Info<< "Selecting an IQNILSRelaxation model for " << dict.dictName()
        << " with initial relaxation factor " << this->initRelax_ << "\n"
        << "\twhile reusing coupling data from " << this->reuse_ 
        << " time steps"
        << endl;
}

template<class Type>
Foam::IQNILSRelaxation<Type>::IQNILSRelaxation
(
    const IQNILSRelaxation<Type>& fR
)
:
    relaxationModel<Type>(fR),
    reuse_(fR.reuse_),
    V_(fR.V_),
    W_(fR.W_),
    T_(fR.T_),
    fldRef_(fR.fldRef_),
    resRef_(fR.resRef_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::IQNILSRelaxation<Type>::~IQNILSRelaxation()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::IQNILSRelaxation<Type>::initialize(const Field<Type> &curFld)
{
    Foam::relaxationModel<Type>::initialize(curFld);
}

template<class Type>
void Foam::IQNILSRelaxation<Type>::relax(Field<Type> &curFld)
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

    updateVW(curFld);

    if (T_.size() > 1)
    {
        Info<< nl
            << "Relaxing field with IQN-ILS procedure" 
            << endl;

        relaxIQNILS(curFld);
    }
    else
    {   
        Info<< nl
            << "Relaxing field with initial relaxation factor: "
            << this->initRelax_ 
            << endl;

        //- Relax with a field relaxation factor
        curFld = this->prevFld_ + this->initRelax_ * this->resFld_;
    }

    //- Store relaxed field as new field
    this->prevFld_ = curFld;

    //- Increment corrector step counter
    this->corr_++;
}

template<class Type>
void Foam::IQNILSRelaxation<Type>::updateVW(Field<Type> &curFld)
{
    if (this->corr_ == 1)
    {
        Info << "Modes before clean-up : " << T_.size();

        while (true)
        {
            if (T_.size())
            {
                if
                (
                    this->runTime_.timeIndex() - reuse_ 
                  > T_[0]
                )
                {
                    for (label i = 0; i < T_.size()-1; i++)
                    {
                        T_[i] = T_[i+1];
                        V_[i] = V_[i+1];
                        W_[i] = W_[i+1];
                    }
                    T_.remove();
                    V_.remove();
                    W_.remove();
                }
                else
                {
                    break;
                }
            }
            else
            {
                break;
            }
        }

        Info<< ", modes after clean-up : "
            << T_.size() << endl;
    }
    else if (this->corr_ == 2)
    {
        //- Set reference in the first coupling iteration
        fldRef_ = curFld;
        resRef_ = this->resFld_;
    }
    else
    {
        //- Reference has been set in the first coupling iteration
        V_.append
        (
            this->resFld_ - resRef_
        );

        W_.append
        (
            curFld - fldRef_
        );

        T_.append
        (
            this->runTime_.timeIndex()
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "IQNILSRelaxationTemplates.C"
#endif

// ************************************************************************* //
