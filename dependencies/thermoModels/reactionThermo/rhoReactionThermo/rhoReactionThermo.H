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
    Foam::rhoReactionThermo

Description
    Foam::rhoReactionThermo

SourceFiles
    rhoReactionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef rhoReactionThermo_H
#define rhoReactionThermo_H

#include "rhoThermo.H"
#include "basicSpecieMixture.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class rhoReactionThermo Declaration
\*---------------------------------------------------------------------------*/

class rhoReactionThermo
:
    public rhoThermo
{

protected:


public:

    //- Runtime type information
    TypeName("rhoReactionThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoReactionThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Constructors

        //- Construct from mesh and phase name
        rhoReactionThermo(const fvMesh&, const word& phaseName);


    // Selectors

        // static autoPtr<rhoReactionThermo> New(const fvMesh& mesh)
        // {
        //     return New(mesh, mesh);
        // }

        //- Standard selection based on fvMesh
        static autoPtr<rhoReactionThermo> New(const fvMesh&, const word& phaseName=word::null);

        // //- Select and check that package contains 'thermoType'
        // static autoPtr<rhoReactionThermo> NewType
        // (
        //     const fvMesh&,
        //     const word& phaseName,
        //     const word& thermoType
        // );


    //- Destructor
    virtual ~rhoReactionThermo();


    // Member functions

        //- Return the composition of the multi-component mixture
        virtual basicSpecieMixture& composition() = 0;

        //- Return the composition of the multi-component mixture
        virtual const basicSpecieMixture& composition() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
