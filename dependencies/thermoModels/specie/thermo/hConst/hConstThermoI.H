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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class equationOfState>
inline Foam::hConstThermo<equationOfState>::hConstThermo
(
    const equationOfState& st,
    const scalar cp,
    const scalar hf
)
:
    equationOfState(st),
    Cp_(cp),
    Hf_(hf)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::hConstThermo<equationOfState>::hConstThermo
(
    const word& name,
    const hConstThermo& ct
)
:
    equationOfState(name, ct),
    Cp_(ct.Cp_),
    Hf_(ct.Hf_)
{}


template<class equationOfState>
inline Foam::autoPtr<Foam::hConstThermo<equationOfState> >
Foam::hConstThermo<equationOfState>::clone() const
{
    return autoPtr<hConstThermo<equationOfState>>::New(*this);
}

template<class equationOfState>
inline Foam::autoPtr<Foam::hConstThermo<equationOfState> >
Foam::hConstThermo<equationOfState>::New(const dictionary& dict)
{
    return autoPtr<hConstThermo<equationOfState>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}

template<class equationOfState>
inline Foam::scalar Foam::hConstThermo<equationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    return Cp_ + equationOfState::Cp(p, T);
}


template<class equationOfState>
inline Foam::scalar Foam::hConstThermo<equationOfState>::Ha
(
    const scalar p, const scalar T
) const
{
    return Hs(p, T) + Hc();
}


template<class equationOfState>
inline Foam::scalar Foam::hConstThermo<equationOfState>::Hs
(
    const scalar p, const scalar T
) const
{
    // Info << "Print Cp_ in hConstThermoI.H function calculated " << Cp(p, T) << endl;
    // Info << "Print Cp_ in hConstThermoI.H function " << this->Cp_ << endl;
    return Cp_*T + equationOfState::H(p, T);
}


template<class equationOfState>
inline Foam::scalar Foam::hConstThermo<equationOfState>::Hc() const
{
    return Hf_;
}


template<class equationOfState>
inline Foam::scalar Foam::hConstThermo<equationOfState>::S
(
    const scalar p, const scalar T
) const
{
    return Cp_*log(T/this->Tstd) + equationOfState::S(p, T);
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline void Foam::hConstThermo<equationOfState>::operator+=
(
    const hConstThermo<equationOfState>& ct
)
{
    // scalar Y1 = this->Y();

    // equationOfState::operator+=(ct);

    // Y1 /= this->Y();
    // scalar Y2 = ct.Y()/this->Y();

    // Cp_ = Y1*Cp_ + Y2*ct.Cp_;
    // Hf_ = Y1*Hf_ + Y2*ct.Hf_;

    scalar Y1 = this->Y();

    equationOfState::operator+=(ct);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        scalar Y2 = ct.Y()/this->Y();

        Cp_ = Y1*Cp_ + Y2*ct.Cp_;
        Hf_ = Y1*Hf_ + Y2*ct.Hf_;
    }
}


// template<class equationOfState>
// inline void Foam::hConstThermo<equationOfState>::operator-=
// (
//     const hConstThermo<equationOfState>& ct
// )
// {
//     scalar Y1 = this->Y();

//     equationOfState::operator-=(ct);

//     Y1 /= this->Y();
//     scalar Y2 = ct.Y()/this->Y();

//     Cp_ = Y1*Cp_ - Y2*ct.Cp_;
//     Hf_ = Y1*Hf_ - Y2*ct.Hf_;
// }


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::hConstThermo<equationOfState> Foam::operator+
(
    const hConstThermo<equationOfState>& ct1,
    const hConstThermo<equationOfState>& ct2
)
{
    // equationOfState eofs
    // (
    //     static_cast<const equationOfState&>(ct1)
    //   + static_cast<const equationOfState&>(ct2)
    // );

    // return hConstThermo<equationOfState>
    // (
    //     eofs,
    //     ct1.Y()/eofs.Y()*ct1.Cp_
    //   + ct2.Y()/eofs.Y()*ct2.Cp_,
    //     ct1.Y()/eofs.Y()*ct1.Hf_
    //   + ct2.Y()/eofs.Y()*ct2.Hf_
    // );

    equationOfState eofs
    (
        static_cast<const equationOfState&>(ct1)
      + static_cast<const equationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < SMALL)
    {
        return hConstThermo<equationOfState>
        (
            eofs,
            ct1.Cp_,
            ct1.Hf_
        );
    }
    else
    {
        return hConstThermo<equationOfState>
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.Cp_
          + ct2.Y()/eofs.Y()*ct2.Cp_,
            ct1.Y()/eofs.Y()*ct1.Hf_
          + ct2.Y()/eofs.Y()*ct2.Hf_
        );
    }
}


// template<class equationOfState>
// inline Foam::hConstThermo<equationOfState> Foam::operator-
// (
//     const hConstThermo<equationOfState>& ct1,
//     const hConstThermo<equationOfState>& ct2
// )
// {
//     equationOfState eofs
//     (
//         static_cast<const equationOfState&>(ct1)
//       - static_cast<const equationOfState&>(ct2)
//     );

//     return hConstThermo<equationOfState>
//     (
//         eofs,
//         ct1.Y()/eofs.Y()*ct1.Cp_
//       - ct2.Y()/eofs.Y()*ct2.Cp_,
//         ct1.Y()/eofs.Y()*ct1.Hf_
//       - ct2.Y()/eofs.Y()*ct2.Hf_
//     );
// }


template<class equationOfState>
inline Foam::hConstThermo<equationOfState> Foam::operator*
(
    const scalar s,
    const hConstThermo<equationOfState>& ct
)
{
    return hConstThermo<equationOfState>
    (
        s*static_cast<const equationOfState&>(ct),
        ct.Cp_,
        ct.Hf_
    );
}


template<class equationOfState>
inline Foam::hConstThermo<equationOfState> Foam::operator==
(
    const hConstThermo<equationOfState>& ct1,
    const hConstThermo<equationOfState>& ct2
)
{
    // return ct2 - ct1;

    equationOfState eofs
    (
        static_cast<const equationOfState&>(ct1)
     == static_cast<const equationOfState&>(ct2)
    );

    return hConstThermo<equationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.Cp_
      - ct1.Y()/eofs.Y()*ct1.Cp_,
        ct2.Y()/eofs.Y()*ct2.Hf_
      - ct1.Y()/eofs.Y()*ct1.Hf_
    );
}


// ************************************************************************* //
