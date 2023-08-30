/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
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

#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Construct from components
inline specie::specie
(
    const word& name,
    const scalar Y,
    const scalar molWeight
)
:
    name_(name),
    Y_(Y),
    molWeight_(molWeight)
{}


// Construct from components without name
inline specie::specie
(
    const scalar Y,
    const scalar molWeight
)
:
    Y_(Y),
    molWeight_(molWeight)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// // Construct as copy
// inline specie::specie(const specie& st)
// :
//     name_(st.name_),
//     Y_(st.Y_),
//     molWeight_(st.molWeight_)
// {}


// Construct as named copy
inline specie::specie(const word& name, const specie& st)
:
    name_(name),
    Y_(st.Y_),
    molWeight_(st.molWeight_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const word& specie::name() const
{
    return name_;
}

//- Molecular weight [kg/kmol]
inline scalar specie::W() const
{
    return molWeight_;
}

//- Molecular volume of specie [m^3/mol]
inline scalar specie::V() const
{
    return molarVolume_;
}
//- Valence of specie (i.e. charge number)
inline scalar specie::z() const
{
    return z_;
}

//- No of moles of this species in mixture
inline scalar specie::Y() const
{
    return Y_;
}

//- Gas constant [J/(kg K)]
inline scalar specie::R() const
{
    return RR/molWeight_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void specie::operator=(const specie& st)
{
    // Info << "Print name of the specie in operator= in specie " << st.name() << endl;  // Is normal that this "doesn't work"
    Y_ = st.Y_;   // It seems to be normal that moles are not 1 here
    molWeight_ = st.molWeight_;
    // Info << "Print Y_ " << Y_ << endl;
    // Info << "Print molWeight_ " << molWeight_ << endl;  // It seems to be normal that molecular weight as well .. WHY?? (Alos ESi-Version)
    // Info << "Print debug" << endl;
}


inline void specie::operator+=(const specie& st)
{   
    // TODO: I don't get this function
    // scalar sumY = max(Y_ + st.Y_, SMALL);

    // molWeight_ =
    //     Y_/sumY*molWeight_
    //   + st.Y_/sumY*st.molWeight_;

    // Y_ = sumY;

    const scalar sumMoles = Y_ + st.Y_;
    if (mag(sumMoles) > SMALL)
    {
        molWeight_ = sumMoles/(Y_/molWeight_ + st.Y_/st.molWeight_);
    }

    Y_ = sumMoles;
}


// inline void specie::operator-=(const specie& st)
// {
//     scalar diffY_ = Y_ - st.Y_;
//     if (mag(diffY_) < SMALL)
//     {
//         diffY_ = SMALL;
//     }

//     molWeight_ =
//         Y_/diffY_*molWeight_
//       - st.Y_/diffY_*st.molWeight_;

//     Y_ = diffY_;
// }


inline void specie::operator*=(const scalar s)
{
    Y_ *= s;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline specie operator+(const specie& st1, const specie& st2)
{
    // scalar sumY = max(st1.Y_ + st2.Y_, SMALL);

    // return specie
    // (
    //     sumY,
    //     st1.Y_/sumY*st1.molWeight_
    //   + st2.Y_/sumY*st2.molWeight_
    // );

    const scalar sumY = max(st1.Y_ + st2.Y_, SMALL);

    if (mag(sumY) > SMALL)
    {
        return specie
        (
            sumY,
            sumY/(st1.Y_/st1.molWeight_ + st2.Y_/st2.molWeight_)
        );
    }
    else
    {
        return st1;
    }
}


// inline specie operator-(const specie& st1, const specie& st2)
// {
//     scalar diffY_ = st1.Y_ - st2.Y_;
//     if (mag(diffY_) < SMALL)
//     {
//         diffY_ = SMALL;
//     }

//     return specie
//     (
//         diffY_,
//         st1.Y_/diffY_*st1.molWeight_
//       - st2.Y_/diffY_*st2.molWeight_
//     );
// }


inline specie operator*(const scalar s, const specie& st)
{
    return specie
    (
        s*st.Y_,
        st.molWeight_
    );
}


inline specie operator==(const specie& st1, const specie& st2)
{
    // return st2 - st1;

    scalar diffMoles = st2.Y_ - st1.Y_;
    if (mag(diffMoles) < SMALL)
    {
        diffMoles = SMALL;
    }

    const scalar diffRW = st2.Y_/st2.molWeight_ - st1.Y_/st1.molWeight_;

    #ifdef __clang__
    // Using intermediate volatile bool to prevent compiler optimising out the
    // if block (above) - CLANG 3.7.1
    volatile const bool valid = (mag(diffRW) > SMALL);
    const scalar molWeight = valid ? diffMoles/diffRW : GREAT;
    #else
    scalar molWeight = GREAT;
    if (mag(diffRW) > SMALL)
    {
         molWeight = diffMoles/diffRW;
    }
    #endif

    return specie(diffMoles, molWeight);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
