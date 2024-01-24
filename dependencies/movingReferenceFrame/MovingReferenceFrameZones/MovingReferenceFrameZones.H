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
    Foam::MovingReferenceFrameZoness

Description
    Abstract base class for MovingReferenceFrameZoness.

SourceFiles
    MovingReferenceFrameZones.C

\*---------------------------------------------------------------------------*/

#ifndef MovingReferenceFrameZones_H
#define MovingReferenceFrameZones_H

#include "dictionary.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedVector.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "fvMatrices.H"

#include "autoPtr.H"
#include "tmp.H"
#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

#include "IOPtrList.H"

#include "translationalMRFZones.H"
#include "rotationalMRFZones.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class MovingReferenceFrameZones Declaration
\*---------------------------------------------------------------------------*/

class MovingReferenceFrameZones
{

protected:

    // Protected data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Pointer list of translational MRF zones
        translationalMRFZones translationalMRFs_;

        //- Pointer list of rotational MRF zones
        rotationalMRFZones rotationalMRFs_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        MovingReferenceFrameZones(const MovingReferenceFrameZones&);

        //- Disallow default bitwise assignment
        void operator=(const MovingReferenceFrameZones&);

public:

    // Constructors

        //- Construct from fvMesh and Istream
        MovingReferenceFrameZones(const fvMesh& mesh);

    //- Destructor

        virtual ~MovingReferenceFrameZones()
        {}

    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        rotationalMRFZones& rotationalMRFs()
        {
            return rotationalMRFs_;
        }

        translationalMRFZones& translationalMRFs()
        {
            return translationalMRFs_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
