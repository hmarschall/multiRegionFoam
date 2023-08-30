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
    Foam::currentTransferInterface.C

SourceFiles
    currentTransferInterface.C

Author
    Constantin Habes (constantin.habes@tu-darmstadt.de, Affiliation A)

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
    
    This file is part of the multiRegionFoam library.

    You may refer to this software as :
    https://doi.org/10.48550/arXiv.2306.01924

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

    This file is part of the multiRegionFoam library.

    You may refer to this software as :
    https://doi.org/10.48550/arXiv.2306.01924

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#ifndef currentTransferInterface_H
#define currentTransferInterface_H

#include "regionInterfaceType.H"
#include "IOReferencer.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionInterfaces
{

/*---------------------------------------------------------------------------*\
                           Class currentTransferInterface Declaration
\*---------------------------------------------------------------------------*/

class currentTransferInterface
:
    public regionInterfaceType
{

    // Private data

        // interface dictionary
        dictionary dict_;

        // Name of fluid phase name
        word fluidPhaseName_;

        // Name of the fluid patch name
        word fluidPatchName_;

        // Sign of the overpotential
        // // - Return for now the magnitude from the fluid regions
        scalar etaSign_;

public:

    TypeName("currentTransferInterface");

    // Constructors

        currentTransferInterface
        (
            const word& type,
            const dictionary& dict,
            const Time& runTime,
            const fvPatch& patchA,
            const fvPatch& patchB
        );

    // Destructor
    virtual ~currentTransferInterface()
    {}

    const dimensionedScalar& K() const;

    const word fluidPhaseName() const
    {
        return fluidPhaseName_;
    }

    const word fluidPatchName() const
    {
        return fluidPatchName_;
    }

    const scalar etaSign() const
    {
        return etaSign_;
    }

    // Correct interface physics
    void correct();

    // Get min time step based on stability criteria
    scalar getMinDeltaT();

    // Log interface specific information
    void info() const
    {};
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionInterfaces
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
