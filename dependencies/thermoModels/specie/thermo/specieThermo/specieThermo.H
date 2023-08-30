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

Class
    Foam::specieThermo

Description
    Basic thermodynamics type based on the use of fitting functions for
    cp, h, s obtained from the template argument type thermo.  All other
    properties are derived from these primitive functions.

SourceFiles
    specieThermoI.H
    specieThermo.C

\*---------------------------------------------------------------------------*/

#ifndef specieThermo_H
#define specieThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// #include "optimisationSwitch.H"
// #include "tolerancesSwitch.H"

namespace Foam
{
namespace species // Encapsulate it from possible mixture values
{

// Forward declaration of friend functions and operators

template<class Thermo, template<class> class Type> class specieThermo;

template<class Thermo, template<class> class Type>
inline specieThermo<Thermo, Type> operator+
(
    const specieThermo<Thermo, Type>&,
    const specieThermo<Thermo, Type>&
);

// template<class Thermo>
// inline specieThermo<Thermo> operator-
// (
//     const specieThermo<Thermo>&,
//     const specieThermo<Thermo>&
// );

template<class Thermo, template<class> class Type>
inline specieThermo<Thermo, Type> operator*
(
    const scalar,
    const specieThermo<Thermo, Type>&
);

template<class Thermo, template<class> class Type>
inline specieThermo<Thermo, Type> operator==
(
    const specieThermo<Thermo, Type>&,
    const specieThermo<Thermo, Type>&
);

template<class Thermo, template<class> class Type>
Ostream& operator<<
(
    Ostream&,
    const specieThermo<Thermo, Type>&
);


/*---------------------------------------------------------------------------*\
                           Class specieThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, template<class> class Type>
class specieThermo
:
    public Thermo,
    public Type<specieThermo<Thermo, Type>>
{
    // Private data

        // //- Convergence tolerance of energy -> temperature inversion functions
        // static const debug::tolerancesSwitch tol_;

        // //- Max temperature jump of energy -> temperature inversion functions
        // static const debug::tolerancesSwitch TJump_;

        // //- Max number of iterations in energy->temperature inversion functions
        // static const debug::optimisationSwitch maxIter_;

        //- Convergence tolerance of energy -> temperature inversion functions
        static const scalar tol_;

        //- Max number of iterations in energy->temperature inversion functions
        static const int maxIter_;



    // Private member functions

        //- Return the temperature corresponding to the value of the
        //  thermodynamic property f, given the function f = F(T) and dF(T)/dT
        // inline scalar T
        // (
        //     scalar f,
        //     scalar T0,
        //     scalar (specieThermo::*F)(const scalar) const,
        //     scalar (specieThermo::*dFdT)(const scalar) const
        // ) const;

        inline scalar T
        (
            scalar f,
            scalar p,
            scalar T0,
            scalar (specieThermo::*F)(const scalar, const scalar) const,
            scalar (specieThermo::*dFdT)(const scalar, const scalar) const,
            scalar (specieThermo::*limit)(const scalar) const
        ) const;


public:

    //- The thermodynamics of the individual species'
    typedef specieThermo<Thermo, Type> thermoType;

    // Constructors

        //- construct from components
        inline specieThermo(const Thermo& sp);

        // //- Construct from Istream
        // specieThermo(Istream&);

        //- Construct from dictionary
        specieThermo(const dictionary& dict);

        //- Construct as named copy
        inline specieThermo(const word& name, const specieThermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            // return thermo::typeName();
            return Thermo::typeName() + ','
            + Type<specieThermo<Thermo, Type>>::typeName();
            // return "species::specieThermo<" + thermo::typeName() + '>';
        }

        //- Name of Enthalpy/Internal energy
        static inline word heName();

        // Fundamental properties
        // (These functions must be provided in derived types)

            // Heat capacity at constant pressure [J/(kmol K)]
            //scalar cp(const scalar) const;

            // Enthalpy [J/kmol]
            //scalar h(const scalar) const;

            // Sensible enthalpy [J/kmol]
            //scalar hs(const scalar) const;

            // Chemical enthalpy [J/kmol]
            //scalar hc(const scalar) const;

            // Entropy [J/(kmol K)]
            //scalar s(const scalar) const;

        // Calculate and return derived properties
        // (These functions need not provided in derived types)

    // Mass and molar specific properties

            //- Heat capacity at constant pressure/volume [J/(kg K)] -> is actually the sensible enthalpy
            inline scalar Cpv(const scalar p, const scalar T) const;

            //- Enthalpy/Internal energy [J/kg] -> is actually the sensible enthalpy
            inline scalar HE(const scalar p, const scalar T) const;

            //- Ratio of heat capacity at constant pressure to that at
            //  constant pressure/volume []
            inline scalar CpByCpv(const scalar p, const scalar T) const;

        //- Heat capacity at constant pressure/volume

            //- Heat capacity at constant pressure [J/(kmol K)]
            inline scalar cp(const scalar p, const scalar T) const;

            // //- Heat capacity at constant volume [J/(kmol K)]
            // inline scalar Cv(const scalar p, const scalar T) const;

            //- Heat capacity at constant volume [J/(kg K)]
            inline scalar cv(const scalar p, const scalar T) const;


        //- Gibbs enthalpy
            //- Gibbs free energy [J/kg]
            inline scalar G(const scalar p, const scalar T) const;

            //- Gibbs free energy [J/kmol]
            inline scalar g(const scalar p, const scalar T) const;

        //- Helmholtz free energy     
            //- Helmholtz free energy [J/kg]
            inline scalar A(const scalar p, const scalar T) const;

            //- Helmholtz free energy [J/kmol]
            inline scalar a(const scalar p, const scalar T) const;

        //- Sensible internal energy
            // //- Sensible internal energy [J/kmol]
            // inline scalar Es(const scalar p, const scalar T) const;

            inline scalar es(const scalar p, const scalar T) const;

        //- Absolute internal energy
            // //- Absolute internal energy [J/kmol]
            // inline scalar Ea(const scalar p, const scalar T) const;

            inline scalar ea(const scalar p, const scalar T) const;


    // Mole specific derived properties

            //- gamma = cp/cv []
            inline scalar gamma(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kmol]
            inline scalar ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kmol]
            inline scalar hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kmol]
            inline scalar hc() const;

            //- Entropy [J/(kmol K)]
            inline scalar s(const scalar p, const scalar T) const;

            //- Enthalpy/Internal energy [J/kmol]
            inline scalar he(const scalar p, const scalar T) const;


    //     // Equilibrium reaction thermodynamics

    //         //- Equilibrium constant [] i.t.o fugacities
    //         //  = PIi(fi/Pstd)^nui
    //         inline scalar K(const scalar p, const scalar T) const;

    //         //- Equilibrium constant [] i.t.o. partial pressures
    //         //  = PIi(pi/Pstd)^nui
    //         //  For low pressures (where the gas mixture is near perfect) Kp = K
    //         inline scalar Kp(const scalar p, const scalar T) const;

    //         //- Equilibrium constant i.t.o. molar concentration
    //         //  = PIi(ci/cstd)^nui
    //         //  For low pressures (where the gas mixture is near perfect)
    //         //  Kc = Kp(pstd/(RR*T))^nu
    //         inline scalar Kc(const scalar p, const scalar T) const;

    //         //- Equilibrium constant [] i.t.o. mole-fractions
    //         //  For low pressures (where the gas mixture is near perfect)
    //         //  Kx = Kp(pstd/p)^nui
    //         inline scalar Kx(const scalar T, const scalar p) const;

    //         //- Equilibrium constant [] i.t.o. number of moles
    //         //  For low pressures (where the gas mixture is near perfect)
    //         //  Kn = Kp(n*pstd/p)^nui where n = number of moles in mixture
    //         inline scalar Kn
    //         (
    //             const scalar T,
    //             const scalar p,
    //             const scalar n
    //         ) const;


        // Energy->temperature  inversion functions

            //- Temperature from Enthalpy given an initial temperature T0
            inline scalar THE
            (
                const scalar H,
                const scalar p,
                const scalar T0
            ) const;

            //- Temperature from sensible Enthalpy given an initial T0
            inline scalar THs
            (
                const scalar Hs,
                const scalar p,
                const scalar T0
            ) const;

            //- Temperature from absolute enthalpy
            //  given an initial temperature T0
            inline scalar THa
            (
                const scalar H,
                const scalar p,
                const scalar T0
            ) const;

            //- Temperature from sensible internal energy
            //  given an initial temperature T0
            inline scalar TEs
            (
                const scalar E,
                const scalar p,
                const scalar T0
            ) const;

            //- Temperature from absolute internal energy
            //  given an initial temperature T0
            inline scalar TEa
            (
                const scalar E,
                const scalar p,
                const scalar T0
            ) const;


    // Member operators

        inline void operator+=(const specieThermo&);
        // inline void operator-=(const specieThermo&);

        inline void operator*=(const scalar);

    // I-O

        //- Write to Ostream
        void write(Ostream& os) const;


    // Friend operators

        friend specieThermo operator+ <Thermo, Type>
        (
            const specieThermo&,
            const specieThermo&
        );

        // friend specieThermo operator- <Thermo, Type>
        // (
        //     const specieThermo&,
        //     const specieThermo&
        // );

        friend specieThermo operator* <Thermo, Type>
        (
            const scalar s,
            const specieThermo&
        );

        friend specieThermo operator== <Thermo, Type>
        (
            const specieThermo&,
            const specieThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo, Type>
        (
            Ostream&,
            const specieThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace species
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "specieThermoI.H"

#ifdef NoRepository
#   include "specieThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
