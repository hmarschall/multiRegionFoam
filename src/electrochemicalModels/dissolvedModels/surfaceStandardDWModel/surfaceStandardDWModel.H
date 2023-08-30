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
    Foam::epRegionModel.C

SourceFiles
    epRegionModel.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceStandardDWModel_H
#define surfaceStandardDWModel_H

#include "surfaceDissolvedModel.H"
#include "HashTable.H"

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace surfaceDissolvedModels
{

/*---------------------------------------------------------------------------*\
                           Class dissolved Declaration
\*---------------------------------------------------------------------------*/

class surfaceStandardDWModel
:
    public surfaceDissolvedModel
{
        //- Dictionary
        const dictionary& dict_;

        //- Water transfer coefficient
        HashTable<scalar> ksi_;

        //- nd
        dimensionedScalar nd_;

        //- rho / membrane equivillent weight (mol/m3)
        dimensionedScalar rhoOnEW_;

        //- Name of current
        word iName_;

        //- Name of temperature
        word TName_;

        //- Relax factor
        scalar relax_;

        //- Correct factor for diffusion coefficient
        scalar corr_;

        //- Functions

        //- diffusion coefficient
        scalar DNaf(const scalar&, const scalar&) const;

        //- Effective water content
        scalar effLam(const scalar&) const;


public:

    //- Runtime type information
    TypeName("surfaceStandard");

    // Constructors

        surfaceStandardDWModel
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~surfaceStandardDWModel();


    // Member Functions
        //- update
        virtual void update(const word& clName);

        // //- update dmdt at the faceCells
        // virtual void updatePatch(const word& patchName);

        //- correct
        virtual void correct();

        //- solve, if necessary
        virtual void solve();

    //- IO

        virtual bool read(const dictionary& dict);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace surfaceStandardDWModels
} // End namespace Foam

#endif

// ************************************************************************* //
