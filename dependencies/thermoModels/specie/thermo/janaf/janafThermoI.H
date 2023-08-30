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

#include "janafThermo.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class equationOfState>
inline Foam::janafThermo<equationOfState>::janafThermo
(
    const equationOfState& st,
    const scalar Tlow,
    const scalar Thigh,
    const scalar Tcommon,
    const typename janafThermo<equationOfState>::coeffArray& highCpCoeffs,
    const typename janafThermo<equationOfState>::coeffArray& lowCpCoeffs,
    const bool convertCoeffs
)
:
    equationOfState(st),
    Tlow_(Tlow),
    Thigh_(Thigh),
    Tcommon_(Tcommon)
{
    if (convertCoeffs)
    {
        for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
        {
            highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel]*this->R();
            lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel]*this->R();
        }
    }
    else
    {
        for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
        {
            highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel];
            lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel];
        }
    }
}


template<class equationOfState>
inline const typename Foam::janafThermo<equationOfState>::coeffArray&
Foam::janafThermo<equationOfState>::coeffs
(
    const scalar T
) const
{
    if (T < Tcommon_)
    {
        return lowCpCoeffs_;
    }
    else
    {
        return highCpCoeffs_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::janafThermo<equationOfState>::janafThermo
(
    const word& name,
    const janafThermo& jt
)
:
    equationOfState(name, jt),
    Tlow_(jt.Tlow_),
    Thigh_(jt.Thigh_),
    Tcommon_(jt.Tcommon_)
{
    for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        highCpCoeffs_[coefLabel] = jt.highCpCoeffs_[coefLabel];
        lowCpCoeffs_[coefLabel] = jt.lowCpCoeffs_[coefLabel];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::scalar Foam::janafThermo<equationOfState>::limit
(
    const scalar T
) const
{
    if (T < Tlow_ || T > Thigh_)
    {
        WarningInFunction
            << "attempt to use janafThermo<equationOfState>"
               " out of temperature range "
            << Tlow_ << " -> " << Thigh_ << ";  T = " << T
            << endl;

        return min(max(T, Tlow_), Thigh_);
    }
    else
    {
        return T;
    }
}


template<class equationOfState>
inline Foam::scalar Foam::janafThermo<equationOfState>::Tlow() const
{
    return Tlow_;
}


template<class equationOfState>
inline Foam::scalar Foam::janafThermo<equationOfState>::Thigh() const
{
    return Thigh_;
}


template<class equationOfState>
inline Foam::scalar Foam::janafThermo<equationOfState>::Tcommon() const
{
    return Tcommon_;
}


template<class equationOfState>
inline const typename Foam::janafThermo<equationOfState>::coeffArray&
Foam::janafThermo<equationOfState>::highCpCoeffs() const
{
    return highCpCoeffs_;
}


template<class equationOfState>
inline const typename Foam::janafThermo<equationOfState>::coeffArray&
Foam::janafThermo<equationOfState>::lowCpCoeffs() const
{
    return lowCpCoeffs_;
}


template<class equationOfState>
inline Foam::scalar Foam::janafThermo<equationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return
        ((((a[4]*T + a[3])*T + a[2])*T + a[1])*T + a[0])
      + equationOfState::Cp(p, T);
}


template<class equationOfState>
inline Foam::scalar Foam::janafThermo<equationOfState>::Ha
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return
    (
        ((((a[4]/5.0*T + a[3]/4.0)*T + a[2]/3.0)*T + a[1]/2.0)*T + a[0])*T
      + a[5]
    ) + equationOfState::H(p, T);
}


template<class equationOfState>
inline Foam::scalar Foam::janafThermo<equationOfState>::Hs
(
    const scalar p,
    const scalar T
) const
{
    return Ha(p, T) - Hc();
}


template<class equationOfState>
inline Foam::scalar Foam::janafThermo<equationOfState>::Hc() const
{
    const coeffArray& a = lowCpCoeffs_;
    const scalar Tstd = Foam::specie::Tstd;
    return
    (
        (
            (((a[4]/5.0*Tstd + a[3]/4.0)*Tstd + a[2]/3.0)*Tstd + a[1]/2.0)*Tstd
          + a[0]
        )*Tstd + a[5]
    );
}


template<class equationOfState>
inline Foam::scalar Foam::janafThermo<equationOfState>::S
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return
    (
        (((a[4]/4.0*T + a[3]/3.0)*T + a[2]/2.0)*T + a[1])*T + a[0]*log(T)
      + a[6]
    ) + equationOfState::S(p, T);
}


// template<class equationOfState>
// inline Foam::scalar Foam::janafThermo<equationOfState>::dGdT
// (
//     const scalar p,
//     const scalar T
// ) const
// {
//     const coeffArray& a = coeffs(T);
//     return -((a[0] + a[5]/T)/T + a[1]/2 + T*(a[2]/3 + T*(a[3]/4 + T*a[4]/5)));
// }


// template<class equationOfState>
// inline Foam::scalar Foam::janafThermo<equationOfState>::dCpdT
// (
//     const scalar p,
//     const scalar T
// ) const
// {
//     const coeffArray& a = coeffs(T);
//     return
//         (((4*a[4]*T + 3*a[3])*T + 2*a[2])*T + a[1]);
// }

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline void Foam::janafThermo<equationOfState>::operator+=
(
    const janafThermo<equationOfState>& jt
)
{
    scalar Y1 = this->Y();

    equationOfState::operator+=(jt);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = jt.Y()/this->Y();

        Tlow_ = max(Tlow_, jt.Tlow_);
        Thigh_ = min(Thigh_, jt.Thigh_);

        if
        (
            janafThermo<equationOfState>::debug
         && notEqual(Tcommon_, jt.Tcommon_)
        )
        {
            FatalErrorInFunction
                << "Tcommon " << Tcommon_ << " for "
                << (this->name().size() ? this->name() : "others")
                << " != " << jt.Tcommon_ << " for "
                << (jt.name().size() ? jt.name() : "others")
                << exit(FatalError);
        }

        for
        (
            label coefLabel=0;
            coefLabel<janafThermo<equationOfState>::nCoeffs_;
            coefLabel++
        )
        {
            highCpCoeffs_[coefLabel] =
                Y1*highCpCoeffs_[coefLabel]
              + Y2*jt.highCpCoeffs_[coefLabel];

            lowCpCoeffs_[coefLabel] =
                Y1*lowCpCoeffs_[coefLabel]
              + Y2*jt.lowCpCoeffs_[coefLabel];
        }
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::janafThermo<equationOfState> Foam::operator+
(
    const janafThermo<equationOfState>& jt1,
    const janafThermo<equationOfState>& jt2
)
{
    equationOfState eofs = jt1;
    eofs += jt2;

    if (mag(eofs.Y()) < SMALL)
    {
        return janafThermo<equationOfState>
        (
            eofs,
            jt1.Tlow_,
            jt1.Thigh_,
            jt1.Tcommon_,
            jt1.highCpCoeffs_,
            jt1.lowCpCoeffs_
        );
    }
    else
    {
        const scalar Y1 = jt1.Y()/eofs.Y();
        const scalar Y2 = jt2.Y()/eofs.Y();

        typename janafThermo<equationOfState>::coeffArray highCpCoeffs;
        typename janafThermo<equationOfState>::coeffArray lowCpCoeffs;

        for
        (
            label coefLabel=0;
            coefLabel<janafThermo<equationOfState>::nCoeffs_;
            coefLabel++
        )
        {
            highCpCoeffs[coefLabel] =
                Y1*jt1.highCpCoeffs_[coefLabel]
              + Y2*jt2.highCpCoeffs_[coefLabel];

            lowCpCoeffs[coefLabel] =
                Y1*jt1.lowCpCoeffs_[coefLabel]
              + Y2*jt2.lowCpCoeffs_[coefLabel];
        }

        if
        (
            janafThermo<equationOfState>::debug
         && notEqual(jt1.Tcommon_, jt2.Tcommon_)
        )
        {
            FatalErrorInFunction
                << "Tcommon " << jt1.Tcommon_ << " for "
                << (jt1.name().size() ? jt1.name() : "others")
                << " != " << jt2.Tcommon_ << " for "
                << (jt2.name().size() ? jt2.name() : "others")
                << exit(FatalError);
        }

        return janafThermo<equationOfState>
        (
            eofs,
            max(jt1.Tlow_, jt2.Tlow_),
            min(jt1.Thigh_, jt2.Thigh_),
            jt1.Tcommon_,
            highCpCoeffs,
            lowCpCoeffs
        );
    }
}


template<class equationOfState>
inline Foam::janafThermo<equationOfState> Foam::operator*
(
    const scalar s,
    const janafThermo<equationOfState>& jt
)
{
    return janafThermo<equationOfState>
    (
        s*static_cast<const equationOfState&>(jt),
        jt.Tlow_,
        jt.Thigh_,
        jt.Tcommon_,
        jt.highCpCoeffs_,
        jt.lowCpCoeffs_
    );
}


template<class equationOfState>
inline Foam::janafThermo<equationOfState> Foam::operator==
(
    const janafThermo<equationOfState>& jt1,
    const janafThermo<equationOfState>& jt2
)
{
    equationOfState eofs
    (
        static_cast<const equationOfState&>(jt1)
     == static_cast<const equationOfState&>(jt2)
    );

    const scalar Y1 = jt2.Y()/eofs.Y();
    const scalar Y2 = jt1.Y()/eofs.Y();

    typename janafThermo<equationOfState>::coeffArray highCpCoeffs;
    typename janafThermo<equationOfState>::coeffArray lowCpCoeffs;

    for
    (
        label coefLabel=0;
        coefLabel<janafThermo<equationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        highCpCoeffs[coefLabel] =
            Y1*jt2.highCpCoeffs_[coefLabel]
          - Y2*jt1.highCpCoeffs_[coefLabel];

        lowCpCoeffs[coefLabel] =
            Y1*jt2.lowCpCoeffs_[coefLabel]
          - Y2*jt1.lowCpCoeffs_[coefLabel];
    }

    if
    (
        janafThermo<equationOfState>::debug
     && notEqual(jt2.Tcommon_, jt1.Tcommon_)
    )
    {
        FatalErrorInFunction
            << "Tcommon " << jt2.Tcommon_ << " for "
            << (jt2.name().size() ? jt2.name() : "others")
            << " != " << jt1.Tcommon_ << " for "
            << (jt1.name().size() ? jt1.name() : "others")
            << exit(FatalError);
    }

    return janafThermo<equationOfState>
    (
        eofs,
        max(jt2.Tlow_, jt1.Tlow_),
        min(jt2.Thigh_, jt1.Thigh_),
        jt2.Tcommon_,
        highCpCoeffs,
        lowCpCoeffs
    );
}


// ************************************************************************* //
