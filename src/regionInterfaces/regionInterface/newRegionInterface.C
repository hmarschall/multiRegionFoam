/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Fork:     foam-extend
    \\  /    A nd           | Version:  4.0                                 
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
    regionInterface

SourceFiles
    regionInterface.C

Author
    Holger Marschall (holger.marschall@tu-darmstadt.de, Affiliation A)

Contact
    Holger Marschall (holger.marschall@tu-darmstadt.de)
    main developer and principal investigator
    TU Darmstadt

Affiliations
    Affiliation A)
    Computational Multiphase Flow
    Research Group Leader
    Department of Mathematics
    Technical University of Darmstadt, Germany

Acknowledgement

Description
    
    This file is part of the multiRegion coupling library.

    You may refer to this software as :
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#include "regionInterface.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<regionInterface> regionInterface::New
(
    const word& type,
    const dictionary& dict,
    const Time& runTime,
    const fvPatch& patchA,
    const fvPatch& patchB
)
{
    IOdictionaryConstructorTable::iterator cstrIter =
        IOdictionaryConstructorTablePtr_->find(type);

    if (cstrIter == IOdictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "interface::New()"
        )   << "Unknown interface type "
            << type << endl << endl
            << "Valid regionInterfaceTypes are : " << endl
            << IOdictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<regionInterface>
        (cstrIter()(type, dict, runTime, patchA, patchB));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
