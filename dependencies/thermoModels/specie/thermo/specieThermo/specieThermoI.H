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

#include "specieThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::species::specieThermo<Thermo, Type>::specieThermo
(
    const Thermo& sp
)
:
    Thermo(sp)
{}


// template<class Thermo, template<class> class Type>
// inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::T
// (
//     scalar f,
//     scalar T0,
//     scalar (specieThermo<Thermo, Type>::*F)(const scalar) const,
//     scalar (specieThermo<Thermo, Type>::*dFdT)(const scalar) const
// ) const
// {

//     Info << "Print Hallo specieThermo " << endl;
//     scalar Test = T0;
//     scalar Tnew = T0;
//     scalar Ttol = T0*tol_();
//     int iter = 0;

//     do
//     {
//         // Limit the temperature jump in a single corrector to TJump_
//         // HJ, 12/Oct/2010

//         Test = Tnew;
//         Tnew = Test
//             - Foam::min(((this->*F)(Test) - f)/(this->*dFdT)(Test), TJump_());

//         if (iter++ > maxIter_())
//         {
//             // Improvements: graceful exit with recovery.  HJ, 11/Oct/2010
//             InfoIn
//             (
//                 "specieThermo<Thermo, Type>::T(scalar f, scalar T0, "
//                 "scalar (specieThermo<Thermo, Type>::*F)(const scalar) const, "
//                 "scalar (specieThermo<Thermo, Type>::*dFdT)(const scalar) const"
//                 ") const"
//             )   << "Maximum number of iterations exceeded.  T0 = "
//                 << T0 << " F = " << (this->*F)(Test)
//                 << " dFdT = " << (this->*dFdT)(Test)
//                 << endl;

//             // Use value where dFdT is calculated using T0.  HJ, 11/Oct/2010
//             Tnew = f/(this->*dFdT)(T0);

//             return Tnew;
//         }

//     } while (mag(Tnew - Test) > Ttol);

//     return Tnew;
// }


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::T
(
    scalar f,
    scalar p,
    scalar T0,
    scalar (specieThermo<Thermo, Type>::*F)(const scalar, const scalar) const,
    scalar (specieThermo<Thermo, Type>::*dFdT)(const scalar, const scalar) const,
    scalar (specieThermo<Thermo, Type>::*limit)(const scalar) const
) const
{
    if (T0 < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: " << T0
            << abort(FatalError);
    }

    scalar Test = T0;
    scalar Tnew = T0;
    scalar Ttol = T0*tol_;
    int    iter = 0;

    do
    {
        Test = Tnew;
        Tnew =
            (this->*limit)
            (Test - ((this->*F)(p, Test) - f)/(this->*dFdT)(p, Test));

        if (iter++ > maxIter_)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter_
                << abort(FatalError);
        }

    } while (mag(Tnew - Test) > Ttol);

    return Tnew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::species::specieThermo<Thermo, Type>::specieThermo
(
    const word& name,
    const specieThermo& st
)
:
    Thermo(name, st)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::word
Foam::species::specieThermo<Thermo, Type>::heName()
{
    return Type<specieThermo<Thermo, Type>>::energyName();
}

template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::specieThermo<Thermo, Type>::Cpv(const scalar p, const scalar T) const
{
    return Type<specieThermo<Thermo, Type>>::Cpv(*this, p, T);
}

template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::specieThermo<Thermo, Type>::CpByCpv
(
    const scalar p,
    const scalar T
) const
{
    return Type<specieThermo<Thermo, Type>>::CpByCpv(*this, p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::specieThermo<Thermo, Type>::HE(const scalar p, const scalar T) const
{
    return Type<specieThermo<Thermo, Type>>::HE(*this, p, T);
}

//- Heat capacity at constant pressure/volume
//  - Cp is calculated in thermo
template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::cp(const scalar p, const scalar T) const
{
    return this->Cp(p, T)*this->W();
}

// template<class Thermo, template<class> class Type>
// inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::Cv(const scalar p, const scalar T) const
// {
//     return this->Cp(p, T) - this->RR*this->W();
// }

template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::cv(const scalar p, const scalar T) const
{
    return this->Cv(p, T)*this->W();
}


//- Gibbs enthalpy
template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::specieThermo<Thermo, Type>::G(const scalar p, const scalar T) const
{
    return this->Ha(p, T) - T*this->S(p, T);
}

template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::specieThermo<Thermo, Type>::g(const scalar p, const scalar T) const
{
    return this->G(p, T)*this->W();
}


//- Helmholtz free energy
template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::specieThermo<Thermo, Type>::A(const scalar p, const scalar T) const
{
    Info << "Try to call for A (Helmholtz free energy): Need to be implemented as follows " << endl; 
    // return this->Ea(p, T) - T*this->S(p, T);
    return this->Ha(p, T) - p/this->rho(p, T);
}

template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::specieThermo<Thermo, Type>::a(const scalar p, const scalar T) const
{
        Info << "Try to call for a (Helmholtz free energy): Need to be implemented as follows " << endl; 
    // return this->Ea(p, T) - T*this->S(p, T);
    return A(p, T)*this->W();
}

// //- Sensible internal energy
// template<class Thermo, template<class> class Type>
// inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::Es(const scalar p, const scalar T) const
// {
//     return this->Hs(p, T) - p/this->rho(p, T);
// }

template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::es(const scalar p, const scalar T) const
{
    return this->Es(p, T)*this->W();
}

// //- Absolute internal energy
// template<class Thermo, template<class> class Type>
// inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::Ea(const scalar p, const scalar T) const
// {
//     return this->Ha(p, T) - p/this->rho(p, T);
// }

template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::ea(const scalar p, const scalar T) const
{
    return this->Ea(p, T)*this->W();
}


// Mole specific derived properties

//- gamma = cp/cv []
template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::gamma(const scalar p, const scalar T) const
{
    // return cp(p, T) / cv(p, T);
    #ifdef __clang__
    volatile const scalar Cp = this->Cp(p, T);
    #else
    const scalar Cp = this->Cp(p, T);
    #endif

    return Cp/(Cp - this->CpMCv(p, T));
}

//- Absolute Enthalpy
template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::ha(const scalar p, const scalar T) const
{
    return this->Ha(p, T)*this->W();
}

template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::hs(const scalar p, const scalar T) const
{
    return this->Hs(p, T)*this->W();
}

template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::hc() const
{
    return this->Hc()*this->W();
}

template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::s(const scalar p, const scalar T) const
{
    return this->S(p, T)*this->W();
}

template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::he(const scalar p, const scalar T) const
{
    return this->HE(p, T)*this->W();
}



// template<class Thermo, template<class> class Type>
// inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::K(const scalar p, const scalar T) const
// {
//     scalar arg = -this->Y()*this->G(this->Pstd, T)/(this->RR*T);

//     if (arg < 600.0)
//     {
//         return ::exp(arg);
//     }
//     else
//     {
//         return VGREAT;
//     }
// }


// template<class Thermo, template<class> class Type>
// inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::Kp(const scalar p, const scalar T) const
// {
//     return K(p, T);
// }


// template<class Thermo, template<class> class Type>
// inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::Kc(const scalar p, const scalar T) const
// {
//     const scalar nm = this->Y()/this->W();

//     if (equal(nm, SMALL))
//     {
//         return Kp(p, T);
//     }
//     else
//     {
//         return Kp(p, T)*pow(this->Pstd/(this->RR*T), nm);
//     }
// }


// template<class Thermo, template<class> class Type>
// inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::Kx
// (
//     const scalar T,
//     const scalar p
// ) const
// {

//     const scalar nm = this->Y()/this->W();

//     if (equal(nm, SMALL))
//     {
//         return Kp(p, T);
//     }
//     else
//     {
//         return Kp(p, T)*pow(this->Pstd()/p, nm);
//     }
// }


// template<class Thermo, template<class> class Type>
// inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::Kn
// (
//     const scalar T,
//     const scalar p,
//     const scalar n
// ) const
// {

//     const scalar nm = this->Y()/this->W();

//     if (equal(nm, SMALL))
//     {
//         return Kp(p, T);
//     }
//     else
//     {
//         return Kp(p, T)*pow(n*this->Pstd()/p, nm);
//     }
// }

// Energy->temperature  inversion functions

template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::THE
(
    const scalar he,
    const scalar p,
    const scalar T0
) const
{
    // const scalar ths = THs(he, p, T0);
    // return ths;
    return Type<specieThermo<Thermo, Type>>::THE(*this, he, p, T0);
}

template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::THs
(
    const scalar hs,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        hs,
        p,
        T0,
        &specieThermo<Thermo, Type>::Hs,
        &specieThermo<Thermo, Type>::Cp,
        &specieThermo<Thermo, Type>::limit
    );
}

template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::THa
(
    const scalar ha,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        ha,
        p,
        T0,
        &specieThermo<Thermo, Type>::Ha,
        &specieThermo<Thermo, Type>::Cp,
        &specieThermo<Thermo, Type>::limit
    );
}

template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::TEs
(
    const scalar es,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        es,
        p,
        T0,
        &specieThermo<Thermo, Type>::Es,
        &specieThermo<Thermo, Type>::Cv,
        &specieThermo<Thermo, Type>::limit
    );
}

template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::specieThermo<Thermo, Type>::TEa
(
    const scalar ea,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        ea,
        p,
        T0,
        &specieThermo<Thermo, Type>::Ea,
        &specieThermo<Thermo, Type>::Cv,
        &specieThermo<Thermo, Type>::limit
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline void Foam::species::specieThermo<Thermo, Type>::operator+=
(
    const specieThermo<Thermo, Type>& st
)
{
    Thermo::operator+=(st);
}

// template<class Thermo, template<class> class Type>
// inline void Foam::species::specieThermo<Thermo, Type>::operator-=
// (
//     const specieThermo<Thermo, Type>& st
// )
// {
//     Thermo::operator-=(st);
// }

template<class Thermo, template<class> class Type>
inline void Foam::species::specieThermo<Thermo, Type>::operator*=(const scalar s)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::species::specieThermo<Thermo, Type> Foam::species::operator+
(
    const specieThermo<Thermo, Type>& st1,
    const specieThermo<Thermo, Type>& st2
)
{
    return specieThermo<Thermo, Type>
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );
}


// template<class Thermo, template<class> class Type>
// inline Foam::species::specieThermo<Thermo, Type> Foam::species::operator-
// (
//     const specieThermo<Thermo, Type>& st1,
//     const specieThermo<Thermo, Type>& st2
// )
// {
//     return specieThermo<Thermo, Type>
//     (
//         static_cast<const Thermo&>(st1) - static_cast<const Thermo&>(st2)
//     );
// }


template<class Thermo, template<class> class Type>
inline Foam::species::specieThermo<Thermo, Type> Foam::species::operator*
(
    const scalar s,
    const specieThermo<Thermo, Type>& st
)
{
    return specieThermo<Thermo, Type>
    (
        s*static_cast<const Thermo&>(st)
    );
}


template<class Thermo, template<class> class Type>
inline Foam::species::specieThermo<Thermo, Type> Foam::species::operator==
(
    const specieThermo<Thermo, Type>& st1,
    const specieThermo<Thermo, Type>& st2
)
{
    // return st2 - st1;
        return specieThermo<Thermo, Type>
    (
        static_cast<const Thermo&>(st1) == static_cast<const Thermo&>(st2)
    );
}


// ************************************************************************* //
