/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
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

Class
    Foam::capillaryInterface.C

SourceFiles
    capillaryInterface.C

\*---------------------------------------------------------------------------*/

#ifndef capillaryInterface_H
#define capillaryInterface_H

#include "regionInterface.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFieldsFwd.H"
#include "IOReferencer.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionInterfaces
{

/*---------------------------------------------------------------------------*\
                           Class capillaryInterface Declaration
\*---------------------------------------------------------------------------*/

class capillaryInterface
:
    public regionInterface
{
    // Private data

        // interface dictionary
        dictionary dict_;

        // Surface tension
        dimensionedScalar sigma0_;

        areaScalarField sigma_;

protected:

    // Private data

        //- surface velocity field
        mutable autoPtr<areaVectorField> UsPtr_;

        //- surface flux field
        mutable autoPtr<edgeScalarField> phisPtr_;

private:

    // Private Member Functions

        //- Make surface velocity
        void makeUs() const;

        //- Make surface flux
        void makePhis() const;

protected:

    // Protected Member Functions

        //- Clear all demand-driven data
        virtual void clearOut() const;

public:

    TypeName("capillaryInterface");

    // Constructors

        capillaryInterface
        (
            const word& type,
            const dictionary& dict,
            const Time& runTime,
            const fvPatch& patchA,
            const fvPatch& patchB
        );

    // Destructor
    virtual ~capillaryInterface()
    {}

    // Public Member Functions

        // Access

        // Correct interface physics - here surface tension
        void correct();

        // Get min time step based on stability criteria
        scalar getMinDeltaT();

        // Edit

        // Update surface velocity
        void updateUs();

        // Update surface flux
        void updatePhis();

        //- Correct surface velocity b.c.
        void correctUsBoundaryConditions();

        //- interfacial transport

        //-- return surface velocity
        Foam::areaVectorField& Us()
        {
            if (UsPtr_.empty())
            {
                makeUs();
            }

            return UsPtr_();
        }

        const Foam::areaVectorField& Us() const
        {
            if (UsPtr_.empty())
            {
                makeUs();
            }

            return UsPtr_();
        }

        //-- return surface flux
        edgeScalarField& Phis()
        {
            if (phisPtr_.empty())
            {
                makePhis();
            }

            return phisPtr_();
        }

        const edgeScalarField& Phis() const
        {
            if (phisPtr_.empty())
            {
               makePhis();
            }

            return phisPtr_();
        }

        //- interfacial forces

        const areaScalarField& sigma() const
        {
            return sigma_;
        }

        tmp<vectorField> surfaceTensionForce() const;

        tmp<vectorField> tangentialSurfaceTensionForce() const;

        tmp<vectorField> normalSurfaceTensionForce() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionInterfaces
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
