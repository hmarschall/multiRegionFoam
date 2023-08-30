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
    Foam::eConstThermo

Group
    grpSpecieThermo

Description
    Constant properties thermodynamics package templated on an equation of
    state

SourceFiles
    eConstThermoI.H
    eConstThermo.C

\*---------------------------------------------------------------------------*/

#ifndef eConstThermo_H
#define eConstThermo_H

#include "specieThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class equationOfState> class eConstThermo;

template<class equationOfState>
inline eConstThermo<equationOfState> operator+
(
    const eConstThermo<equationOfState>&,
    const eConstThermo<equationOfState>&
);

template<class equationOfState>
inline eConstThermo<equationOfState> operator*
(
    const scalar,
    const eConstThermo<equationOfState>&
);

template<class equationOfState>
inline eConstThermo<equationOfState> operator==
(
    const eConstThermo<equationOfState>&,
    const eConstThermo<equationOfState>&
);

template<class equationOfState>
Ostream& operator<<
(
    Ostream&,
    const eConstThermo<equationOfState>&
);


/*---------------------------------------------------------------------------*\
                        Class eConstThermo Declaration
\*---------------------------------------------------------------------------*/

template<class equationOfState>
class eConstThermo
:
    public equationOfState
{
    // Private data

        //- Heat capacity at constant volume
        //  Note: input in [J/(kg K)], but internally uses [J/(kmol K)]
        scalar Cv_;

        //- Heat of formation
        //  Note: input in [J/kg], but internally uses [J/kmol]
        scalar Hf_;


    // Private Member Functions

        //- Construct from components
        inline eConstThermo
        (
            const equationOfState& st,
            const scalar cv,
            const scalar hf
        );


public:

    // Constructors

        //- Construct from dictionary
        eConstThermo(const dictionary& dict);

        //- Construct as named copy
        inline eConstThermo(const word&, const eConstThermo&);

        //- Construct and return a clone
        inline autoPtr<eConstThermo> clone() const;

        // Selector from dictionary
        inline static autoPtr<eConstThermo> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "eConst<" + equationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Fundamental properties

            //- Heat capacity at constant volume [J/(kg K)]
            inline scalar Cv(const scalar p, const scalar T) const;

            //- Sensible internal energy [J/kg]
            inline scalar Es(const scalar p, const scalar T) const;

            //- Absolute internal energy [J/kg]
            inline scalar Ea(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            //- Entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

        // EtoHthermo
            #include "EtoHthermo.H"
            

        // // Derivative term used for Jacobian

        //     //- Derivative of Gibbs free energy w.r.t. temperature
        //     inline scalar dGdT(const scalar p, const scalar T) const;

        //     //- Temperature derivative of heat capacity at constant pressure
        //     inline scalar dCpdT(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const eConstThermo&);


    // Friend operators

        friend eConstThermo operator+ <equationOfState>
        (
            const eConstThermo&,
            const eConstThermo&
        );

        friend eConstThermo operator* <equationOfState>
        (
            const scalar,
            const eConstThermo&
        );

        friend eConstThermo operator== <equationOfState>
        (
            const eConstThermo&,
            const eConstThermo&
        );


    // IOstream Operators

        friend Ostream& operator<< <equationOfState>
        (
            Ostream&,
            const eConstThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "eConstThermoI.H"

#ifdef NoRepository
    #include "eConstThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
