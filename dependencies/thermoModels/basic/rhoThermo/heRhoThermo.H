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
    Foam::heRhoThermo

Description
    Foam::heRhoThermo

SourceFiles
    heRhoThermo.C

\*---------------------------------------------------------------------------*/

#ifndef heRhoThermo_H
#define heRhoThermo_H

// #include "rhoReactionThermo.H"
#include "rhoThermo.H"
#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class heRhoThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicThermo, class MixtureType>
class heRhoThermo
:
    public heThermo<BasicThermo, MixtureType>
{
    // // Private member functions

    //     // void calculate();

        //- Calculate the thermo variables
        void calculate
        (
            const volScalarField& p,
            volScalarField& T,
            volScalarField& he,
            volScalarField& psi,
            volScalarField& rho,
            volScalarField& mu,
            volScalarField& alpha,
            const bool doOldTimes
        );

        //- Construct as copy (not implemented)
        heRhoThermo(const heRhoThermo<BasicThermo, MixtureType>&);


public:

    //- Runtime type information
    TypeName("heRhoThermo");


    // Constructors

    //- Construct from mesh
    heRhoThermo(const fvMesh&, const word&);        // obj is the mesh (see hREactionThermo.C)


    //- Destructor
    virtual ~heRhoThermo();


    // Member functions

        //- Update properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "heRhoThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
