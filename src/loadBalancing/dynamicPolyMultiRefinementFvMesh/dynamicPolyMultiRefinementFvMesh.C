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

\*---------------------------------------------------------------------------*/

#include "dynamicPolyMultiRefinementFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"
#include "refinementSelection.H"
#include "enhancedPrismatic2DRefinement.H"
#include "enhancedPolyhedralRefinement.H"
#include "compositeRefinementSelection.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dynamicPolyMultiRefinementFvMesh, 0);

addToRunTimeSelectionTable
(
    topoChangerFvMesh,
    dynamicPolyMultiRefinementFvMesh,
    IOobject
);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicPolyMultiRefinementFvMesh::dynamicPolyMultiRefinementFvMesh
(
    const IOobject& io,
    const word subDictName
)
:
    topoChangerFvMesh(io),
    refinementDict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).subDict(subDictName + "Coeffs")
    ),
    refinements_(),
    refinementInterval_(readLabel(refinementDict_.lookup("refineInterval"))),
    unrefinementInterval_(readLabel(refinementDict_.lookup("unrefineInterval")))
{
    // Check whether we read polyMeshModifiers from
    // constant/polyMesh/meshModifiers file in the base class
    if (!topoChanger_.empty())
    {
        // Already initialized, warn the user that we'll neglect it
        WarningInFunction
           << "Using controls from constant/dynamicMeshDict instead of"
           << " constant/polyMesh/meshModifiers."
           << nl
           << "To supress this warning, delete meshModifiers file."
           << endl;

       // Clear the list
       topoChanger_.clear();
    }

    // Only one topo changer engine
    topoChanger_.setSize(1);

    // Get number of valid geometric dimensions
    const label nGeometricDirs = this->nGeometricD();

    switch(nGeometricDirs)
    {
        case 3:
            // Add the polyhedralRefinement engine for
            // 3D isotropic refinement
            Info<< "3D case detected. "
                << "Adding enhancedPolyhedralRefinement topology modifier" << endl;
            topoChanger_.set
            (
                0,
                new enhancedPolyhedralRefinement
                (
                    "enhancedPolyhedralRefinement",
                    refinementDict_,
                    0,
                    topoChanger_
                )
            );
            break;

        case 2:
            // Add the enhancedPrismatic2DRefinement engine for
            // 2D isotropic refinement
            Info<< "2D case detected. "
                << "Adding enhancedPrismatic2DRefinement topology modifier" << endl;
            topoChanger_.set
            (
                0,
                new enhancedPrismatic2DRefinement
                (
                    "enhancedPrismatic2DRefinement",
                    refinementDict_,
                    0,
                    topoChanger_
                )
            );
            break;

        case 1:
            FatalErrorInFunction
                << "1D case detected. No valid refinement strategy is"
                <<  " available for 1D cases."
                << abort(FatalError);
            break;

        default:
            FatalErrorInFunction
                << "Invalid number of geometric meshes detected: "
                << nGeometricDirs
                << nl << "It appears that this mesh is neither 1D, 2D or 3D."
                << abort(FatalError);

    }

    // Write mesh and modifiers
    topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
    topoChanger_.write();
    write();

    // Initialize refinement selection algorithm after modifiers
    // refinement selection as entries
    Info<< "Reading refinements list." << endl;
    const PtrList<entry> refinsInfo ( refinementDict_.lookup("refinements") );
    if (refinsInfo.empty()) {
        FatalIOErrorInFunction(refinementDict_)
            << "At least one refinement request is expected."
            << exit(FatalIOError);
    }
    // Reshape refinements list
    refinements_.setSize(refinsInfo.size());

    // Construct refinement list
    forAll(refinements_, ri)
    {
        // Select a requested refinement
        const entry& refinInfo = refinsInfo[ri];
        // Require that the refinement entry is a valid dict
        if(!refinInfo.isDict())
        {
            FatalIOErrorInFunction(refinementDict_)
                << "Entry " << refinInfo << " in refinements section is not a"
                << " valid dictionary." << exit(FatalIOError);
        }
        dictionary refDict = refinInfo.dict() | refinementDict_;
        refinements_[ri].first() = refinementDetails(time(), io, *this, refDict);
        refinements_[ri].second() = refinementSelection::New(*this, refDict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicPolyMultiRefinementFvMesh::~dynamicPolyMultiRefinementFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicPolyMultiRefinementFvMesh::firstUpdate() const
{
    forAll(refinements_, ri)
    {
        // DISCUSSION: If any refinement thinks this is not the first update
        if (refinements_[ri].first().curTimeIndex() > time().timeIndex())
        {
            return false;
        }
    }
    return true;
}

bool Foam::dynamicPolyMultiRefinementFvMesh::singleRefUpdate()
{
    const label timeID = time().timeIndex();
    bool changed = false;

    // Get reference to base class refinement polyMeshModifier
    enhancedRefinement& refModifier = refCast<enhancedRefinement>(topoChanger_[0]);

    // 1. apply all refinements

    auto& refinData = refinements_[0].first();
    auto& refinSelectionPtr = refinements_[0].second();

    // Re-read the data from dictionary for on-the-fly changes
    refinData.readDict();

    if 
    (
        (timeID % refinementInterval_) == 0
        && refinData.curTimeIndex() < timeID
    )
    {
        // Update current time index to skip multiple topo changes per single
        // time step
        refinData.setCurTimeIndex( time().timeIndex() );
        refinData.setRefinementSettings(refModifier);

        // Create empty list for refinement candidates
        labelList refCandidates;

        refCandidates.transfer
        (
            refinSelectionPtr->refinementCellCandidates()()
        );
        Info<< "Selected " << refCandidates.size()
            << " refinement candidates."
            << endl;

        // Set cells to refine. Note: refinement needs to make sure that face
        // and point consistent refinement is performed
        refModifier.setCellsToRefine(refCandidates);

        labelList unrefCandidates;
        if (timeID > 2 && (timeID % unrefinementInterval_) == 0)
        {
            unrefCandidates.transfer
            (
                refinSelectionPtr->unrefinementPointCandidates()()
            );
            refModifier.setSplitPointsToUnrefine(unrefCandidates);
        }

        //// Activate the polyhedral refinement engine if there are some cells to
        //// refine or there are some split points to unrefine around
        bool enableTopoChange = !refCandidates.empty() || !unrefCandidates.empty();

        //// Note: must enable topo change for all processors since face and point
        //// consistent refinement must be ensured across coupled patches
        reduce(enableTopoChange, orOp<bool>());

        if (enableTopoChange)
        {
            refModifier.enable();
        }
        else
        {
            refModifier.disable();
        }

        // Perform refinement and unrefinement in one go
        if (!polyTopoChanger::debug) Pout.stdStream().setstate(std::ios_base::failbit);
        autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();
        if (!polyTopoChanger::debug) Pout.stdStream().clear();

        // Output cell balance if the topo change has been performed
        const label nOldCells =
            returnReduce(topoChangeMap->nOldCells(), sumOp<label>());
        const label sizeCellMap =
            returnReduce(topoChangeMap->cellMap().size(), sumOp<label>());

        // If the size of cell map is different than zero, we actually performed
        // some topo changes
        if (sizeCellMap)
        {
            Info<< "Successfully performed polyhedral refinement/unrefinement. "
                << "Changed from " << nOldCells << " to " << sizeCellMap
                << " cells." << endl;
        }

        changed = topoChangeMap->morphing();
    }
    else
    {
        // Update current time index to skip multiple topo change checks
        // per time step
        refinData.setCurTimeIndex( time().timeIndex() );
    }
    return changed;
}

bool Foam::dynamicPolyMultiRefinementFvMesh::multiRefUpdate()
{
    const label timeID = time().timeIndex();
    bool changed = false;
    List<bool> refining(refinements_.size(), false);
    // If multiple refinements are requested:
    // 0. unrefine around over-refined points
    // 1. apply refinements in succession

    // Get reference to base class refinement polyMeshModifier
    enhancedRefinement& refModifier = refCast<enhancedRefinement>(topoChanger_[0]);

    // 0. Over refined cells
    //labelListList refCandidates(refinements_.size());
    if (timeID > 2)
    {
        boolList overRefinedCell(nCells(), false);
        labelList refCandidate(nCells(), -1);

        forAll(refinements_,  ri)
        {
            // Re-read the data from dictionary for on-the-fly changes
            auto& refinData = refinements_[ri].first();
            auto& refinSelectionPtr = refinements_[ri].second();
            refinData.readDict();

            // Reset refinement settings
            refinData.setRefinementSettings(refModifier);

            // Store refinement candidates
            labelList refCandidates;
            refCandidates.transfer
            (
                refinSelectionPtr->refinementCellCandidates()()
            );

            refModifier.setCellsToRefine(refCandidates);
            refCandidates = refModifier.cellsToRefine();

            forAll(refCandidates, rci) {
                refCandidate[refCandidates[rci]] = ri;
            }
        }

        labelList unrefPoints;
        const labelList& cellLevel = refModifier.cellLevel();
        forAll(refinements_,  ri)
        {
            // if a cell is over-refined for one refinement request, and is not a candidate for any other, unrefine it
            // Note: the cell is a left-over from previous refinements, and we assume it's safe to unrefine it
            auto& refinData = refinements_[ri].first();
            forAll(cellLevel, ci) {
                if (cellLevel[ci] >= refinData.maxRefinementLevel() && refCandidate[ci] == -1)
                {
                    unrefPoints.append(cellPoints()[ci]);
                }
            }
        }

        refModifier.setCellsToRefine(labelList());
        refModifier.setSplitPointsToUnrefine(unrefPoints);

        // Activate the polyhedral refinement engine if there are some split points to unrefine around
        bool enableTopoChange = !unrefPoints.empty();

        // Note: must enable topo change for all processors since face and point
        // consistent refinement must be ensured across coupled patches
        reduce(enableTopoChange, orOp<bool>());

        if (enableTopoChange)
        {
            refModifier.enable();
        }
        else
        {
            refModifier.disable();
        }

        // Perform unrefinement
        if (!polyTopoChanger::debug) Pout.stdStream().setstate(std::ios_base::failbit);
        autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();
        if (!polyTopoChanger::debug) Pout.stdStream().clear();

        // Output cell balance if the topo change has been performed
        const label nOldCells =
            returnReduce(topoChangeMap->nOldCells(), sumOp<label>());
        const label sizeCellMap =
            returnReduce(topoChangeMap->cellMap().size(), sumOp<label>());

        // If the size of cell map is different than zero, we actually performed
        // some topo changes
        if (sizeCellMap)
        {
            Info<< "Successfully performed polyhedral unrefinement. "
                << "Changed from " << nOldCells << " to " << sizeCellMap
                << " cells." << endl;
        }

        changed = topoChangeMap->morphing() ? true : changed;
    }

    // 1. apply all refinements

    forAll(refinements_, ri)
    {
        auto& refinData = refinements_[ri].first();
        auto& refinSelectionPtr = refinements_[ri].second();

        // Re-read the data from dictionary for on-the-fly changes
        refinData.readDict();

        if (refinData.curTimeIndex() < timeID)
        {
            // Update current time index to skip multiple topo changes per single
            // time step
            refinData.setCurTimeIndex( time().timeIndex() );
            refinData.setRefinementSettings(refModifier);

            // Create empty list for refinement candidates
            labelList refCandidates;

            refCandidates.transfer
            (
                refinSelectionPtr->refinementCellCandidates()()
            );
            Info<< "Selected " << refCandidates.size()
                << " refinement candidates."
                << endl;

            // Set cells to refine. Note: refinement needs to make sure that face
            // and point consistent refinement is performed
            refModifier.setCellsToRefine(refCandidates);

            //// Activate the polyhedral refinement engine if there are some cells to
            //// refine or there are some split points to unrefine around
            bool enableTopoChange = !refCandidates.empty() ; //|| !finalUnrefCandidates.empty();

            //// Note: must enable topo change for all processors since face and point
            //// consistent refinement must be ensured across coupled patches
            reduce(enableTopoChange, orOp<bool>());

            if (enableTopoChange)
            {
                refModifier.enable();
            }
            else
            {
                refModifier.disable();
            }

            // Perform refinement and unrefinement in one go
            if (!polyTopoChanger::debug) Pout.stdStream().setstate(std::ios_base::failbit);
            autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();
            if (!polyTopoChanger::debug) Pout.stdStream().clear();

            // Output cell balance if the topo change has been performed
            const label nOldCells =
                returnReduce(topoChangeMap->nOldCells(), sumOp<label>());
            const label sizeCellMap =
                returnReduce(topoChangeMap->cellMap().size(), sumOp<label>());

            // If the size of cell map is different than zero, we actually performed
            // some topo changes
            if (sizeCellMap)
            {
                Info<< "Successfully performed polyhedral refinement. "
                    << "Changed from " << nOldCells << " to " << sizeCellMap
                    << " cells." << endl;
            }

            refining[ri] = topoChangeMap->morphing();
            
        }
        else
        {
            // Update current time index to skip multiple topo change checks
            // per time step
            refinData.setCurTimeIndex( time().timeIndex() );
        }
    }
    forAll(refining, ri)
    {
        changed = refining[ri] ? refining[ri] : changed;
    }
    return changed;
}


bool Foam::dynamicPolyMultiRefinementFvMesh::unrefUpdate()
{
    // Refines only points which all refinement selections agree on unrefining
    // Get reference to base class refinement polyMeshModifier
    enhancedRefinement& refModifier = refCast<enhancedRefinement>(topoChanger_[0]);

    // Accumulate which refinement selection suggests to unrefine a point
    // unrefPoint[pi] == refinements_.size() means point pi needs to be unrefined
    labelList unrefPoint(nPoints(), 0);

    forAll(refinements_,  ri)
    {
        // Re-read the data from dictionary for on-the-fly changes
        auto& refinData = refinements_[ri].first();
        auto& refinSelectionPtr = refinements_[ri].second();
        refinData.readDict();

        // Reset refinement settings
        refinData.setRefinementSettings(refModifier);

        // Store unrefinement candidates
        labelList unrefCandidates;
        unrefCandidates.transfer
        (
            refinSelectionPtr->unrefinementPointCandidates()()
        );

        forAll(unrefCandidates, upi) {
            unrefPoint[unrefCandidates[upi]] += 1;
        }
    }

    labelList unrefPoints;
    forAll(unrefPoint,  pi)
    {
        // if a point is marked as an unrefinement candidate, add it
        if (unrefPoint[pi] == refinements_.size()) unrefPoints.append(pi);
    }

    refModifier.setCellsToRefine(labelList());
    refModifier.setSplitPointsToUnrefine(unrefPoints);

    // Activate the polyhedral refinement engine if there are some split points to unrefine around
    bool enableTopoChange = !unrefPoints.empty();

    // Note: must enable topo change for all processors since face and point
    // consistent refinement must be ensured across coupled patches
    reduce(enableTopoChange, orOp<bool>());

    if (enableTopoChange)
    {
        refModifier.enable();
    }
    else
    {
        refModifier.disable();
    }

    // Perform unrefinement
    if (!polyTopoChanger::debug) Pout.stdStream().setstate(std::ios_base::failbit);
    autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();
    if (!polyTopoChanger::debug) Pout.stdStream().clear();

    // Output cell balance if the topo change has been performed
    const label nOldCells =
        returnReduce(topoChangeMap->nOldCells(), sumOp<label>());
    const label sizeCellMap =
        returnReduce(topoChangeMap->cellMap().size(), sumOp<label>());

    // If the size of cell map is different than zero, we actually performed
    // some topo changes
    if (sizeCellMap)
    {
        Info<< "Successfully performed polyhedral unrefinement. "
            << "Changed from " << nOldCells << " to " << sizeCellMap
            << " cells." << endl;
    }

    return topoChangeMap->morphing() ? true : false;
}


bool Foam::dynamicPolyMultiRefinementFvMesh::update()
{
    // Get time index
    const label timeID = time().timeIndex();
    bool changed = false;

    // which ops should we perform
    bool doUnrefinement = (timeID % unrefinementInterval_ == 0);
    bool doRefinement = (timeID % refinementInterval_ == 0);

    if (timeID > 2 && (!doRefinement && !doUnrefinement))
    {
        // No AMR op is requested, return immediately
        return false;
    }

    // If only one refinement is requested, fall back to single-refinement
    if (refinements_.size() == 1) return singleRefUpdate();

    if (timeID > 2 && doUnrefinement) {
        // Unrefinement requested
        bool unrefChanged = unrefUpdate();
        changed = unrefChanged ? unrefChanged : changed;
    }

    if (doRefinement) {
        // Some refinements are requested
        bool multiChanged = multiRefUpdate();
        changed = multiChanged ? multiChanged : changed;
    }

    return changed;
}


// ************************************************************************* //
