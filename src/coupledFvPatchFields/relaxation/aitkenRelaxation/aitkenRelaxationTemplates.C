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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::aitkenRelaxation<Type>::updateAitkenFactor(Field<Type> &curFld)
{
    notImplemented
    (
        "aitkenRelaxationTemplates.C\n"
        "void aitkenRelaxation<Type>::updateAitkenFactor\n"
        "(\n"
        "    Field<Type> &curFld\n"
        ")\n"
        "not implemented"
    );
}

template<>
inline void Foam::aitkenRelaxation<vector>::updateAitkenFactor(Field<vector> &curFld)
{
    if (this->corr_ == 1)
    {
        aitkenRelax_ = 
            (aitkenRelax_ < this->initRelax_) 
            ? this->initRelax_ : aitkenRelax_;

        return;
    }

    aitkenRelax_ = 
        -aitkenRelax_
       *(
            gSum
            (
                this->prevResFld_
              & (this->resFld_ - this->prevResFld_)
            )
           /(
                gSum
                (
                    (
                        this->resFld_
                      - this->prevResFld_
                    )
                  & (
                        this->resFld_
                      - this->prevResFld_
                    )
                )
            )
        );

    aitkenRelax_ = mag(aitkenRelax_);


    //- CH: Could also use a maximum relaxation factor insead of 1.0
    //  as it is done in solids4foam AitkenCouplingInterface.C
    if (aitkenRelax_ > 1.0)
    {
        Info<< nl
            << "Aitken relaxation factor greater than 1" << nl
            << "Clipping Aitken relaxation factor to 1"
            << endl;

        aitkenRelax_ = 1.0;
    }
}


template<>
inline void Foam::aitkenRelaxation<scalar>::updateAitkenFactor(Field<scalar> &curFld)
{
    if (this->corr_ == 1)
    {
        aitkenRelax_ = 
            (aitkenRelax_ < this->initRelax_) 
            ? this->initRelax_ : aitkenRelax_;

        return;
    }

    aitkenRelax_ = 
        -aitkenRelax_
       *(
            gSum
            (
                this->prevResFld_
              * (this->resFld_ - this->prevResFld_)
            )
           /(
                gSum
                (
                    (
                        this->resFld_
                      - this->prevResFld_
                    )
                  * (
                        this->resFld_
                      - this->prevResFld_
                    )
                )
            )
        );


    aitkenRelax_ = mag(aitkenRelax_);


    //- CH: Could also use a maximum relaxation factor insead of 1.0
    //  as it is done in solids4foam AitkenCouplingInterface.C
    if (aitkenRelax_ > 1.0)
    {
        Info<< nl
            << "Aitken relaxation factor greater than 1" << nl
            << "Clipping Aitken relaxation factor to 1"
            << endl;

        aitkenRelax_ = 1.0;
    }
}

// ************************************************************************* //
