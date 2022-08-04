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

Author
    Daniel Deising, TU Darmstadt. All rights reserved.
    Milad Bagheri, TU Darmstadt. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "regionsRefinement.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "volPointInterpolation.H"
#include "pointVolInterpolation.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(regionsRefinement, 0);
addToRunTimeSelectionTable
(
    refinementSelection,
    regionsRefinement,
    dictionary
);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionsRefinement::regionsRefinement
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    refinementSelection(mesh, dict),
    mesh_(mesh),
    refinedRegions_
    (PtrList<entry>
    (coeffDict().lookup("regions"))
    ),
    cellPointCellSmoothing_
    (
        coeffDict().lookupOrDefault<Switch>("cellPointCellSmoothing", false)
    )
{
}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::regionsRefinement::~regionsRefinement()
{}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

Foam::Xfer<Foam::labelList>
Foam::regionsRefinement::refinementCellCandidates() const
{
    labelList refinementCandidates(mesh().nCells(),0);

        forAll(refinedRegions_, regionI)
        {
            const entry& region = refinedRegions_[regionI];
            autoPtr<topoSetSource> source =
            topoSetSource::New(region.keyword(), mesh_, region.dict());
        
            cellSet selectedCellSet
                (
                    mesh_,
                    "cellSet",
                    mesh_.nCells()/10+1 //Estimation
                );
            source->applyToSet
                (
                    topoSetSource::NEW,
                    selectedCellSet
                );
        
            const labelList cells = selectedCellSet.toc();
            label minLevel = readLabel(region.dict().lookup("minLevel"));
                forAll(cells, i)
                {
                    const label& cellI = cells[i];
                    refinementCandidates[cellI] = 
                    max(refinementCandidates[cellI], minLevel);
                    refinementCandidates.append(cellI);
                }
        }
        
    Info<< "Selection algorithm " << type() << " selected "
        << returnReduce(refinementCandidates.size(), sumOp<label>())
        << " cells as refinement candidates."
        << endl;

    // Return the list in the Xfer container to prevent copying
    return refinementCandidates.xfer();
}


Foam::Xfer<Foam::labelList>
Foam::regionsRefinement::unrefinementPointCandidates() const
{
    // Create the list for unrefinement candidates
    labelList unrefinementCandidates(0);


    // Print out some information
    Info<< "Selection algorithm " << type() << " selected "
        << returnReduce(unrefinementCandidates.size(), sumOp<label>())
        << " points as unrefinement candidates."
        << endl;

    // Return the list in the Xfer container to prevent copying
    return unrefinementCandidates.xfer();
}

// ************************************************************************* //
