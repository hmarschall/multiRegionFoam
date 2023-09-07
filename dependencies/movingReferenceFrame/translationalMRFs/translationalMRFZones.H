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
    Foam::translationalMRFZones

Description
    Container class for a set of translationalMRFZones with the 
    translationalMRFZone member functions  implemented to loop over 
    the functions for each translationalMRFZone.

SourceFiles
    translationalMRFZones.C

\*---------------------------------------------------------------------------*/

#ifndef translationalMRFZones_H
#define translationalMRFZones_H

#include "translationalMRFZone.H"
#include "IOPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class translationalMRFZones Declaration
\*---------------------------------------------------------------------------*/

class translationalMRFZones
:
    public IOPtrList<translationalMRFZone>
{
    // Private data

        //- Mesh reference
        const fvMesh& mesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        translationalMRFZones(const translationalMRFZones&);

        //- Disallow default bitwise assignment
        void operator=(const translationalMRFZones&);


public:

    // Constructors

        //- Construct from fvMesh
        translationalMRFZones(const fvMesh& mesh);


    // Member Functions

        // Correct MRF properties
        void correctMRF();

        //- Add the frame acceleration contribution to the momentum equation
        void addFrameAcceleration
        (
            fvVectorMatrix& UEqn,
            const volScalarField& rho
        );

        void addFrameAcceleration
        (
            fvVectorMatrix& UEqn
        );

        //- Correct the boundary velocity for the rotation of the MRF region
        void correctBoundaryVelocity
        (
            volVectorField& U,
            surfaceScalarField& phi
        );

        //- Write and restart MRF parameters
        void writeRestart();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
