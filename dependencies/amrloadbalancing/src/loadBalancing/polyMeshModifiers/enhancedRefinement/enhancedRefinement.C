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
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.
    Hrvoje Jasak, Wikki Ltd.
    Mohammed Elwardi Fadeli

Notes
    Generalisation of hexRef8 for polyhedral cells and refactorisation into mesh
    modifier engine.

\*---------------------------------------------------------------------------*/

#include "enhancedRefinement.H"
#include "tmp.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(enhancedRefinement, 0);

    // Note: do not add to run-time selection table since this is abstract base
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::enhancedRefinement::enhancedRefinement
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyTopoChanger& mme
)
:
    refinement(name, dict, index, mme)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::enhancedRefinement::~enhancedRefinement()
{}

Foam::label Foam::enhancedRefinement::faceConsistentRefinement
(
    boolList& cellsToRefine
) const
{
    // Count number of cells that will be added
    label nAddCells = 0;

    // Get necessary mesh data
    const label nFaces = mesh_.nFaces();
    const label nInternalFaces = mesh_.nInternalFaces();

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    // Loop through internal faces and check consistency
    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        // Get owner and neighbour labels
        const label& own = owner[faceI];
        const label& nei = neighbour[faceI];

        // Get owner and neighbour cell levels
        // Note: If the cell is marked for refinement, the level is current
        // level + 1, otherwise it is equal to the current level
        // BugFix; look back instead of looking up; fixes the unwanted refinement
        const label ownLevel =
            cellsToRefine[own] ? cellLevel_[own] : cellLevel_[own]-1;
        const label neiLevel =
            cellsToRefine[nei] ? cellLevel_[nei] : cellLevel_[nei]-1;

        if (ownLevel > (neiLevel + 1))
        {
            // Owner level is higher than neighbour level + 1, neighbour must be
            // marked for refinement
            cellsToRefine[nei] = true;
            ++nAddCells;
        }
        else if (neiLevel > (ownLevel + 1))
        {
            // Neighbour level is higher than owner level + 1, owner must be
            // marked for refinement
            cellsToRefine[own] = true;
            ++nAddCells;
        }
    }

    // Create owner level for boundary faces to prepare for swapping on coupled
    // boundaries
    labelList ownLevel(nFaces - nInternalFaces);
    forAll (ownLevel, i)
    {
        // Get owner of the face and update owner cell levels
        const label& own = owner[i + nInternalFaces];
        ownLevel[i] =
            cellsToRefine[own] ? cellLevel_[own] : cellLevel_[own]-1;
    }

    // Swap boundary face lists (coupled boundary update)
    syncTools::swapBoundaryFaceList(mesh_, ownLevel, false);

    // Note: now the ownLevel list actually contains the neighbouring level
    // (from the other side), use alias (reference) for clarity from now on
    const labelList& neiLevel = ownLevel;

    // Loop through boundary faces
    forAll (neiLevel, i)
    {
        // Get owner of the face and owner level
        const label& own = owner[i + nInternalFaces];
        const label curOwnLevel =
            cellsToRefine[own] ? cellLevel_[own]+1 : cellLevel_[own];

        // Note: we are using more stringent 1:1 consistency across coupled
        // boundaries in order to simplify handling of edge based consistency
        // checks for parallel runs
        // Bugfix related to PLB: Check whether owner is already marked for
        // refinement. Will allow 2:1 consistency across certain processor faces
        // where we have a new processor boundary. VV, 23/Jan/2019.
        if
        (
            (neiLevel[i] > curOwnLevel)
         && !cellsToRefine[own]
        )
        {
            // Neighbour level is higher than owner level, owner must be
            // marked for refinement
            cellsToRefine[own] = true;
            ++nAddCells;
        }

        // Note: other possibility (that owner level is higher than neighbour
        // level) is taken into account on the other side automatically
    }

    // Return number of added cells
    return nAddCells;
}

Foam::label Foam::enhancedRefinement::faceConsistentUnrefinement
(
    boolList& cellsToUnrefine
) const
{
    // Count number of removed cells from unrefinement
    label nRemCells = 0;

    // Get necessary mesh data
    const label nFaces = mesh_.nFaces();
    const label nInternalFaces = mesh_.nInternalFaces();

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    // Loop through internal faces and check consistency
    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        // Get owner and neighbour labels
        const label& own = owner[faceI];
        const label& nei = neighbour[faceI];

        // Get owner and neighbour cell levels
        // Note: If the cell is marked for unrefinement, the level is current
        // level - 1, otherwise it is equal to the current level
        // BugFix; look back instead of looking up; fixes the unwanted refinement
        const label ownLevel =
            cellsToUnrefine[own] ? cellLevel_[own]-1 : cellLevel_[own];
        const label neiLevel =
            cellsToUnrefine[nei] ? cellLevel_[nei]-1 : cellLevel_[nei];

        if (ownLevel < (neiLevel - 1))
        {
            // Owner level is smaller than neighbour level - 1, we must not
            // unrefine owner

            // Check whether the cell has not been marked for unrefinement
            if (!cellsToUnrefine[own])
            {
                FatalErrorInFunction
                    << "Cell not marked for unrefinement, indicating a"
                    << " previous unnoticed problem with unrefinement."
                    << nl
                    << "Owner: " << own << ", neighbour: " << nei
                    << nl
                    << "Owner level: " << ownLevel
                    << ", neighbour level: " << neiLevel << nl
                    << "This is probably because the refinement and "
                    << "unrefinement regions are very close." << nl
                    << "Try increasing nUnrefinementBufferLayers. "
                    << abort(FatalError);
            }
            else
            {
                cellsToUnrefine[own] = false;
                ++nRemCells;
            }
        }
        else if (neiLevel < (ownLevel - 1))
        {
            // Neighbour level is smaller than owner level - 1, we must not
            // unrefine neighbour

            // Check whether the cell has not been marked for unrefinement
            if (!cellsToUnrefine[nei])
            {
                FatalErrorInFunction
                    << "Cell not marked for unrefinement, indicating a"
                    << " previous unnoticed problem with unrefinement."
                    << nl
                    << "Owner: " << own << ", neighbour: " << nei
                    << nl
                    << "Owner level: " << ownLevel
                    << ", neighbour level: " << neiLevel << nl
                    << "This is probably because the refinement and "
                    << "unrefinement regions are very close." << nl
                    << "Try increasing nUnrefinementBufferLayers. "
                    << abort(FatalError);
            }
            else
            {
                cellsToUnrefine[nei] = false;
                ++nRemCells;
            }
        }
    }

    // Create owner level for boundary faces to prepare for swapping on coupled
    // boundaries
    labelList ownLevel(nFaces - nInternalFaces);
    forAll (ownLevel, i)
    {
        // Get owner of the face and update owner cell levels
        const label& own = owner[i + nInternalFaces];
        ownLevel[i] =
            cellsToUnrefine[own] ? cellLevel_[own]-1 : cellLevel_[own];
    }

    // Swap boundary face lists (coupled boundary update)
    syncTools::swapBoundaryFaceList(mesh_, ownLevel, false);

    // Note: now the ownLevel list actually contains the neighbouring level
    // (from the other side), use alias (reference) for clarity from now on
    const labelList& neiLevel = ownLevel;

    // Loop through boundary faces
    forAll (neiLevel, i)
    {
        // Get owner of the face and owner level
        const label& own = owner[i + nInternalFaces];
        const label curOwnLevel =
            cellsToUnrefine[own] ? cellLevel_[own] : cellLevel_[own]+1;

        // Note: we are using more stringent 1:1 consistency across coupled
        // boundaries in order to simplify handling of edge based consistency
        // checks for parallel runs
        if (curOwnLevel < neiLevel[i])
        {
            // Owner level is smaller than neighbour level, we must not
            // unrefine owner

            // Check whether the cell has not been marked for unrefinement
            // Note: this "redundancy" check should not be performed if we are
            // running with dynamic load balancing. If an ordinary face with
            // standard consistency (2:1) becomes a processor face with more
            // stringent consistency (1:1), the refinement still remains valid,
            // even though the 1:1 consistency is not achieved for this time
            // step. Doing further refinement will make sure that this does not
            // exceed at least 2:1 consistency and therefore 2:1 edge
            // consistency as well. Instead of issuing a FatalError, issue a
            // Warning and wrap it into debug
            if (!cellsToUnrefine[own])
            {
                if (debug)
                {
                    WarningInFunction
                        << "Boundary cell not marked for unrefinement,"
                        << " indicating a previous unnoticed problem with"
                        << " unrefinement."
                        << nl
                        << "Owner: " << own
                        << nl
                        << "Owner level: " << curOwnLevel
                        << ", neighbour level: " << neiLevel[i] << nl
                        << "This is probably because the refinement and "
                        << "unrefinement regions are very close." << nl
                        << "Try increasing nUnrefinementBufferLayers. "
                        << nl
                        << "Another possibility is that you are running "
                        << "with Dynamic Load Balancing, in which case "
                        << "this should be fine."
                        << endl;
                }
            }
            else
            {
                cellsToUnrefine[own] = false;
                ++nRemCells;
            }
        }

        // Note: other possibility (that neighbour level is smaller than owner
        // level) is taken into account on the other side automatically
    }

    // Return number of local cells removed from unrefinement
    return nRemCells;
}

// ************************************************************************* //
