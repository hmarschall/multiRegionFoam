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

#include "hPolynomialThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class equationOfState, int PolySize>
inline Foam::hPolynomialThermo<equationOfState, PolySize>::hPolynomialThermo
(
    const equationOfState& pt,
    const scalar Hf,
    const scalar Sf,
    const Polynomial<PolySize>& CpCoeffs,
    const typename Polynomial<PolySize>::intPolyType& hCoeffs,
    const Polynomial<PolySize>& sCoeffs
)
:
    equationOfState(pt),
    Hf_(Hf),
    Sf_(Sf),
    CpCoeffs_(CpCoeffs),
    hCoeffs_(hCoeffs),
    sCoeffs_(sCoeffs)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class equationOfState, int PolySize>
inline Foam::hPolynomialThermo<equationOfState, PolySize>::hPolynomialThermo
(
    const word& name,
    const hPolynomialThermo& pt
)
:
    equationOfState(name, pt),
    Hf_(pt.Hf_),
    Sf_(pt.Sf_),
    CpCoeffs_(pt.CpCoeffs_),
    hCoeffs_(pt.hCoeffs_),
    sCoeffs_(pt.sCoeffs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class equationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<equationOfState, PolySize>::limit
(
    const scalar T
) const
{
    return T;
}


template<class equationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<equationOfState, PolySize>::Cp
(
    const scalar p, const scalar T
) const
{
    return CpCoeffs_.value(T) + equationOfState::Cp(p, T);
}


template<class equationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<equationOfState, PolySize>::Ha
(
    const scalar p, const scalar T
) const
{
    return hCoeffs_.value(T) + equationOfState::H(p, T);
}


template<class equationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<equationOfState, PolySize>::Hs
(
    const scalar p, const scalar T
) const
{
    return Ha(p, T) - Hc();
}


template<class equationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<equationOfState, PolySize>::Hc()
const
{
    return Hf_;
}


template<class equationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<equationOfState, PolySize>::S
(
    const scalar p,
    const scalar T
) const
{
    return sCoeffs_.value(T) + equationOfState::S(p, T);
}


// template<class equationOfState, int PolySize>
// inline Foam::scalar Foam::hPolynomialThermo<equationOfState, PolySize>::dGdT
// (
//     const scalar p,
//     const scalar T
// ) const
// {
//     return
//     (
//         hCoeffs_.derivative(T)
//       - T*sCoeffs_.derivative(T)
//       - sCoeffs_.value(T)
//     );
// }


// template<class equationOfState, int PolySize>
// inline Foam::scalar Foam::hPolynomialThermo<equationOfState, PolySize>::dCpdT
// (
//     const scalar p,
//     const scalar T
// ) const
// {
//     return
//     (
//         CpCoeffs_.derivative(T)
//     );
// }

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class equationOfState, int PolySize>
inline void Foam::hPolynomialThermo<equationOfState, PolySize>::operator=
(
    const hPolynomialThermo<equationOfState, PolySize>& pt
)
{
    equationOfState::operator=(pt);

    Hf_ = pt.Hf_;
    Sf_ = pt.Sf_;
    CpCoeffs_ = pt.CpCoeffs_;
    hCoeffs_ = pt.hCoeffs_;
    sCoeffs_ = pt.sCoeffs_;
}


template<class equationOfState, int PolySize>
inline void Foam::hPolynomialThermo<equationOfState, PolySize>::operator+=
(
    const hPolynomialThermo<equationOfState, PolySize>& pt
)
{
    scalar Y = this->Y();

    equationOfState::operator+=(pt);

    if (mag(this->Y()) > SMALL)
    {
        Y /= this->Y();
        const scalar Y2 = pt.Y()/this->Y();

        Hf_ = Y*Hf_ + Y2*pt.Hf_;
        Sf_ = Y*Sf_ + Y2*pt.Sf_;
        CpCoeffs_ = Y*CpCoeffs_ + Y2*pt.CpCoeffs_;
        hCoeffs_ = Y*hCoeffs_ + Y2*pt.hCoeffs_;
        sCoeffs_ = Y*sCoeffs_ + Y2*pt.sCoeffs_;
    }
}


template<class equationOfState, int PolySize>
inline void Foam::hPolynomialThermo<equationOfState, PolySize>::operator*=
(
    const scalar s
)
{
    equationOfState::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class equationOfState, int PolySize>
inline Foam::hPolynomialThermo<equationOfState, PolySize> Foam::operator+
(
    const hPolynomialThermo<equationOfState, PolySize>& pt1,
    const hPolynomialThermo<equationOfState, PolySize>& pt2
)
{
    equationOfState eofs = pt1;
    eofs += pt2;

    if (mag(eofs.Y()) < SMALL)
    {
        return hPolynomialThermo<equationOfState, PolySize>
        (
            eofs,
            pt1.Hf_,
            pt1.Sf_,
            pt1.CpCoeffs_,
            pt1.hCoeffs_,
            pt1.sCoeffs_
        );
    }
    {
        const scalar Y = pt1.Y()/eofs.Y();
        const scalar Y2 = pt2.Y()/eofs.Y();

        return hPolynomialThermo<equationOfState, PolySize>
        (
            eofs,
            Y*pt1.Hf_ + Y2*pt2.Hf_,
            Y*pt1.Sf_ + Y2*pt2.Sf_,
            Y*pt1.CpCoeffs_ + Y2*pt2.CpCoeffs_,
            Y*pt1.hCoeffs_ + Y2*pt2.hCoeffs_,
            Y*pt1.sCoeffs_ + Y2*pt2.sCoeffs_
        );
    }
}


template<class equationOfState, int PolySize>
inline Foam::hPolynomialThermo<equationOfState, PolySize> Foam::operator*
(
    const scalar s,
    const hPolynomialThermo<equationOfState, PolySize>& pt
)
{
    return hPolynomialThermo<equationOfState, PolySize>
    (
        s*static_cast<const equationOfState&>(pt),
        pt.Hf_,
        pt.Sf_,
        pt.CpCoeffs_,
        pt.hCoeffs_,
        pt.sCoeffs_
    );
}


template<class equationOfState, int PolySize>
inline Foam::hPolynomialThermo<equationOfState, PolySize> Foam::operator==
(
    const hPolynomialThermo<equationOfState, PolySize>& pt1,
    const hPolynomialThermo<equationOfState, PolySize>& pt2
)
{
    equationOfState eofs
    (
        static_cast<const equationOfState&>(pt1)
     == static_cast<const equationOfState&>(pt2)
    );

    const scalar Y = pt1.Y()/eofs.Y();
    const scalar Y2 = pt2.Y()/eofs.Y();

    return hPolynomialThermo<equationOfState, PolySize>
    (
        eofs,
        Y2*pt2.Hf_       - Y*pt1.Hf_,
        Y2*pt2.Sf_       - Y*pt1.Sf_,
        Y2*pt2.CpCoeffs_ - Y*pt1.CpCoeffs_,
        Y2*pt2.hCoeffs_  - Y*pt1.hCoeffs_,
        Y2*pt2.sCoeffs_  - Y*pt1.sCoeffs_
    );
}


// ************************************************************************* //
