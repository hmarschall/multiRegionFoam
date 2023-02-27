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

#include "relaxationModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class Type>
Foam::relaxationModel<Type>::relaxationModel
(
    const Time& runTime,
    const dictionary& dict
)
:
    runTime_(runTime),
    curTime_(runTime.value()),
    corr_(0),
    prevFld_(),
    resFld_(),
    prevResFld_(),
    initRelax_(dict.lookupOrDefault<scalar>("relax",1.0))
{}

template<class Type>
Foam::relaxationModel<Type>::relaxationModel
(
    const relaxationModel<Type>& rM
)
:
    runTime_(rM.runTime_),
    curTime_(rM.curTime_),
    corr_(rM.corr_),
    prevFld_(rM.prevFld_),
    resFld_(rM.resFld_),
    prevResFld_(rM.prevResFld_),
    initRelax_(rM.initRelax_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::relaxationModel<Type>::~relaxationModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::relaxationModel<Type>::initialize(const Field<Type> &initFld)
{
    prevFld_ = initFld;
}

template<class Type>
void Foam::relaxationModel<Type>::updateResiual(const Field<Type> &curFld)
{
    //- Store old residual
    prevResFld_ = resFld_;

    //- Compute new residual with curent and previous field
    resFld_ = curFld - prevFld_;
}

template<class Type>
void Foam::relaxationModel<Type>::write(Ostream& os) const
{
    os.writeKeyword("relaxType") << type() << token::END_STATEMENT << nl;
    os.writeKeyword("relax") << initRelax_ << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "newRelaxationModel.C"

// ************************************************************************* //
