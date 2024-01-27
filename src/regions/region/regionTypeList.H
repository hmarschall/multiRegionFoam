/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Fork:     foam-extend
    \\  /    A nd           | Version:  4.1                                 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::regionTypeList

SourceFiles
    regionTypeList.C

Author
    Holger Marschall (holger.marschall@tu-darmstadt.de, Affiliation A)
    Constantin Habes (constantin.habes@tu-darmstadt.de, Affiliation A)
    Heba Alkafri (heba.alkafri@tu-darmstadt.de, Affiliation A)

Contact
    Holger Marschall (holger.marschall@tu-darmstadt.de)
    main developer and principal investigator
    TU Darmstadt

Affiliations
    Affiliation A)
    Computational Multiphase Flow
    Department of Mathematics
    Technical University of Darmstadt, Germany

Acknowledgement
    Funded by
    Hessian Ministry of Higher Education, Research, Science and the Arts,
    National High Performance Computing Center for Computational Engineering
    Science (NHR4CES)

Description
    List container for regions

    This file is part of the multiRegionFoam library.

    You may refer to this software as :
    https://doi.org/10.48550/arXiv.2306.01924

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#ifndef regionTypeList_H
#define regionTypeList_H

#include "dynamicFvMesh.H"
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
        bool usesPIMPLE(const bool warn = false) const;

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

        //- Post solve actions after PIMPLE solution
        void postSolvePIMPLE();

        //- Clear regions
        void clear();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
