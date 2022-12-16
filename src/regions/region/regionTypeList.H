/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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
    Foam::regionTypeList

Description
    List container for regions

SourceFiles
    regionTypeList.C

\*---------------------------------------------------------------------------*/

#ifndef regionTypeList_H
#define regionTypeList_H

#include "dynamicFvMesh.H"
#include "mergePolyMesh.H"
#include "dictionary.H"
#include "hashedWordList.H"
#include "regionProperties.H"
#include "regionType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class regionTypeList Declaration
\*---------------------------------------------------------------------------*/

class regionTypeList
:
    public PtrList<regionType>
{
private:

    // Private Data

        //- dictionary
        const IOdictionary dict_;

        //- runTime
        const Time& runTime_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionTypeList(const regionTypeList&);

        //- Disallow default bitwise assignment
        void operator=(const regionTypeList&);

protected:

    // Protected data

        regionProperties region_;

        //- list of mesh names
        wordList meshNames_;

public:

    //- Constructor
    regionTypeList(const Time& runTime);

    //- Destructor
    virtual ~regionTypeList()
    {}

    // Member Functions

        //- Return active PIMPLE status
        bool active(const bool active = false) const;

        //- Reset the source list
        void reset(const regionProperties& region);

        //- Update the meshes and correct model properties
        void preSolve();

        //- Get minimum time step
        scalar getMinDeltaT();

        //- Solve region-inherent physics
        void solveRegion();

        //- Solve pressure-velocity system using PIMPLE
        void solvePIMPLE();

        //- Correct for mesh motion after PIMPLE loop
        void meshMotionCorrector();

        //- Set equations for coupled solution
        void setCoupledEqns();

        //- Update fields after solution
        void postSolve();

        //- Clear regions
        void clear();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
