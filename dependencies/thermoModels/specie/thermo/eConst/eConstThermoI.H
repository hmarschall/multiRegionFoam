/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class equationOfState>
inline Foam::eConstThermo<equationOfState>::eConstThermo
(
    const equationOfState& st,
    const scalar cv,
    const scalar hf
)
:
    equationOfState(st),
    Cv_(cv),
    Hf_(hf)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::eConstThermo<equationOfState>::eConstThermo
(
    const word& name,
    const eConstThermo<equationOfState>& ct
)
:
    equationOfState(name, ct),
    Cv_(ct.Cv_),
    Hf_(ct.Hf_)
{}


template<class equationOfState>
inline Foam::autoPtr<Foam::eConstThermo<equationOfState>>
Foam::eConstThermo<equationOfState>::clone() const
{
    return autoPtr<eConstThermo<equationOfState>>::New(*this);
}


template<class equationOfState>
inline Foam::autoPtr<Foam::eConstThermo<equationOfState>>
Foam::eConstThermo<equationOfState>::New(const dictionary& dict)
{
    return autoPtr<eConstThermo<equationOfState>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::scalar Foam::eConstThermo<equationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class equationOfState>
inline Foam::scalar Foam::eConstThermo<equationOfState>::Cv
(
    const scalar p,
    const scalar T
) const
{
    return Cv_ + equationOfState::Cv(p, T);
}


template<class equationOfState>
inline Foam::scalar Foam::eConstThermo<equationOfState>::Es
(
    const scalar p,
    const scalar T
) const
{
    return Cv_*T + equationOfState::E(p, T);
}


template<class equationOfState>
inline Foam::scalar Foam::eConstThermo<equationOfState>::Ea
(
    const scalar p,
    const scalar T
) const
{
    return Es(p, T) + Hc();
}


template<class equationOfState>
inline Foam::scalar Foam::eConstThermo<equationOfState>::Hc() const
{
    return Hf_;
}


template<class equationOfState>
inline Foam::scalar Foam::eConstThermo<equationOfState>::S
(
    const scalar p,
    const scalar T
) const
{
    return Cv(p, T)*log(T/this->Tstd) + equationOfState::S(p, T);
}


// template<class equationOfState>
// inline Foam::scalar Foam::eConstThermo<equationOfState>::dGdT
// (
//     const scalar p,
//     const scalar T
// ) const
// {
//     return 0;
// }


// template<class equationOfState>
// inline Foam::scalar Foam::eConstThermo<equationOfState>::dCvdT
// (
//     const scalar p,
//     const scalar T
// ) const
// {
//     return 0;
// }

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline void Foam::eConstThermo<equationOfState>::operator+=
(
    const eConstThermo<equationOfState>& ct
)
{
    scalar Y1 = this->Y();

    equationOfState::operator+=(ct);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        scalar Y2 = ct.Y()/this->Y();

        Cv_ = Y1*Cv_ + Y2*ct.Cv_;
        Hf_ = Y1*Hf_ + Y2*ct.Hf_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::eConstThermo<equationOfState> Foam::operator+
(
    const eConstThermo<equationOfState>& ct1,
    const eConstThermo<equationOfState>& ct2
)
{
    equationOfState eofs
    (
        static_cast<const equationOfState&>(ct1)
      + static_cast<const equationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < SMALL)
    {
        return eConstThermo<equationOfState>
        (
            eofs,
            ct1.Cv_,
            ct1.Hf_
        );
    }
    else
    {
        return eConstThermo<equationOfState>
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.Cv_
          + ct2.Y()/eofs.Y()*ct2.Cv_,
            ct1.Y()/eofs.Y()*ct1.Hf_
          + ct2.Y()/eofs.Y()*ct2.Hf_
        );
    }
}


template<class equationOfState>
inline Foam::eConstThermo<equationOfState> Foam::operator*
(
    const scalar s,
    const eConstThermo<equationOfState>& ct
)
{
    return eConstThermo<equationOfState>
    (
        s*static_cast<const equationOfState&>(ct),
        ct.Cv_,
        ct.Hf_
    );
}


template<class equationOfState>
inline Foam::eConstThermo<equationOfState> Foam::operator==
(
    const eConstThermo<equationOfState>& ct1,
    const eConstThermo<equationOfState>& ct2
)
{
    equationOfState eofs
    (
        static_cast<const equationOfState&>(ct1)
     == static_cast<const equationOfState&>(ct2)
    );

    return eConstThermo<equationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.Cv_
      - ct1.Y()/eofs.Y()*ct1.Cv_,
        ct2.Y()/eofs.Y()*ct2.Hf_
      - ct1.Y()/eofs.Y()*ct1.Hf_
    );
}


// ************************************************************************* //
