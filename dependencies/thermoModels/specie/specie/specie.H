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
    Foam::specie

Description
    Base class of the thermophysical property types.

SourceFiles
    specieI.H
    specie.C

\*---------------------------------------------------------------------------*/

#ifndef specie_H
#define specie_H

#include "word.H"
#include "scalar.H"
#include "debug.H"
#include "constantsSwitch.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class specie;

inline specie operator+(const specie&, const specie&);
inline specie operator*(const scalar, const specie&);
inline specie operator==(const specie&, const specie&);

Ostream& operator<<(Ostream&, const specie&);


/*---------------------------------------------------------------------------*\
                           Class specie Declaration
\*---------------------------------------------------------------------------*/

class specie
{
    // Private data

        //- Name of specie
        word name_;

        //- Mass fraction of this specie in mixture
        scalar Y_;

        //- Molecular weight of specie [kg/kmol]
        scalar molWeight_;

        //- Molecular volume of specie [m^3/mol]
        scalar molarVolume_;

        //- Valence of specie (i.e. charge number)
        scalar z_;


public:

    //- Runtime type information
    ClassName("specie");

    // Public constants

        // Thermodynamic constants

            //- Universal gas constant [J/(mol K)]
            static const scalar RR;

            //- Standard pressure [Pa]
            static const scalar Pstd;

            //- Standard temperature [K]
            static const scalar Tstd;


    // Constructors

        //- Construct from components without name
        inline specie
        (
            const scalar Y,
            const scalar molWeight
        );

        //- Construct from components with name
        inline specie
        (
            const word& name,
            const scalar Y,
            const scalar molWeight
        );

        // //- Construct as copy
        // inline specie(const specie&);

        //- Construct as named copy
        inline specie(const word& name, const specie&);

        //- Construct from dictionary
        specie(const dictionary& dict);


    // Member Functions

        // Access

            //- Name
            inline const word& name() const;

            //- Molecular weight [kg/kmol]
            inline scalar W() const;

            //- Molar volume [m^3/mol]
            inline scalar V() const;

            //- Valence number
            inline scalar z() const;

            //- No of moles of this species in mixture
            inline scalar Y() const;

            //- Gas constant [J/(kg K)]
            inline scalar R() const;

    // IO

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        // This is automatically called when '=' is used with
        // between two specie objects
        inline void operator=(const specie&);
        
        inline void operator+=(const specie&);
        // inline void operator-=(const specie&);
        inline void operator*=(const scalar);


    // Friend operators

        inline friend specie operator+(const specie&, const specie&);
        // inline friend specie operator-(const specie&, const specie&);
        inline friend specie operator*(const scalar, const specie&);
        inline friend specie operator==(const specie&, const specie&);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const specie&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "specieI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
