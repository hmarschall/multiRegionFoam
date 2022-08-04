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
    Holger Marschall, TU Darmstadt.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "constrainedFieldBoundsRefinement.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "volPointInterpolation.H"
#include "pointVolInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(constrainedFieldBoundsRefinement, 0);
addToRunTimeSelectionTable
(
    refinementSelection,
    constrainedFieldBoundsRefinement,
    dictionary
);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constrainedFieldBoundsRefinement::constrainedFieldBoundsRefinement
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    refinementSelection(mesh, dict),
    fieldName_(coeffDict().lookup("fieldName")),
    lowerBound_(readScalar(coeffDict().lookup("lowerBound"))),
    upperBound_(readScalar(coeffDict().lookup("upperBound"))),
    lowerConstrainBound_(readScalar(coeffDict().lookup("lowerConstrainBound"))),
    upperConstrainBound_(readScalar(coeffDict().lookup("upperConstrainBound"))),
    cellPointCellSmoothing_
    (
        coeffDict().lookupOrDefault<Switch>("cellPointCellSmoothing", false)
    )
{}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::constrainedFieldBoundsRefinement::~constrainedFieldBoundsRefinement()
{}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

Foam::Xfer<Foam::labelList>
Foam::constrainedFieldBoundsRefinement::refinementCellCandidates() const
{
    // Get the field
    const volScalarField& vField =
        mesh().lookupObject<volScalarField>(fieldName_);

    // Create temporary for the field (sharing the reference to volume field)
    tmp<volScalarField> tvf(vField);

    // Use cell-point-cell interpolation to smooth out the field
    if (cellPointCellSmoothing_)
    {
        // Create volume to point interpolation object
        const volPointInterpolation& vpi = volPointInterpolation::New(mesh());

        // Interpolate from cell centres to points
        pointScalarField pField(vpi.interpolate(vField));

        // Create point to volume interpolation object
        const pointMesh& pMesh = pointMesh::New(mesh());
        const pointVolInterpolation pvi(pMesh, mesh());

        // Interpolate from points back to cell centres
        tmp<volScalarField> tInterpolatedVf = pvi.interpolate(pField);

        // Assign to the temporary object
        tvf = tInterpolatedVf;
    }

    // Get const reference to (possibly smoothed volume field)
    const volScalarField& vf = tvf();

    // Create storage for collection of cells.
    dynamicLabelList refinementCandidates(Foam::max(100, mesh().nCells()));
    dynamicLabelList excludedRefCandidates(Foam::max(100, mesh().nCells()));

    // Loop through internal field and collect cells to refine
    const scalarField& vfIn = vf.internalField();

    forAll (vfIn, cellI)
    {
        // Get current cell value
        const scalar& cellValue = vfIn[cellI];

        if (cellValue > 0.5)
        {
            // Cell value is within the bounds, append cell for potential
            // refinement
            refinementCandidates.append(cellI);
        }

//        if
//        (
//            cellValue < upperConstrainBound_ 
////         && cellValue > lowerConstrainBound_
//        )
//        {
//            excludedRefCandidates.append(cellI);
//        }
//        else
//        {
//            // Cell value is not within the intervall 
//            // where refinement is constrained
//            if (cellValue > lowerBound_ && cellValue < upperBound_)
//            {
//                // Cell value is within the bounds, append cell for potential
//                // refinement
//                refinementCandidates.append(cellI);
//            }
//        }
    }

    // Print out some information
    Info<< "Selection algorithm " << type() << " selected "
        << returnReduce(refinementCandidates.size(), sumOp<label>())
        << " cells as refinement candidates." 
        << endl;

    // Return the list in the Xfer container to prevent copying
    return refinementCandidates.xfer();
}


Foam::Xfer<Foam::labelList>
Foam::constrainedFieldBoundsRefinement::unrefinementPointCandidates() const
{
    // Get mesh data
    const labelListList& meshCellPoints = mesh().cellPoints();

    // Create mark-up field for points that are found on cells to be unrefined
    boolList pointsOnCellsToRefine(mesh().nPoints(), false);

    // Get the field
    const volScalarField& vField =
        mesh().lookupObject<volScalarField>(fieldName_);

    const scalarField& vFieldIn = vField.internalField();

//    // Create volume to point interpolation object
//    const volPointInterpolation& vpi = volPointInterpolation::New(mesh());

//    // Interpolate the volume field from cell centres to points and get the
//    // internal point field
//    pointScalarField pField(vpi.interpolate(vField));
//    const scalarField& pFieldIn = pField.internalField();

    // Loop through all cells and count number of points on cells to unrefine
    label nPointsUnref = 0;

    // Get current cell level
    const labelIOField& curCellLevel =
        mesh().lookupObject<labelIOField>("cellLevel");

    forAll (vFieldIn, cellI)
    {
        const scalar& cvf = vFieldIn[cellI];
        const scalar& curcl = curCellLevel[cellI];

        if
        (
           (cvf < 0.5)
        )
        {
            // Mark all of points of the cell to be unrefined.
            const labelList& curCellPoints = meshCellPoints[cellI];

            forAll (curCellPoints, i)
            {
                // Get marker for this point
                bool& ptOnRefCell = pointsOnCellsToRefine[curCellPoints[i]];

                if (!ptOnRefCell)
                {
                    // This points has not been marked yet, mark it and
                    // increment the counter for protected points
                    ptOnRefCell = true;
                    ++nPointsUnref;
                }
            }
        }
    }

    // Create the list for unrefinement candidates
    labelList unrefinementCandidates(nPointsUnref);
    label n = 0;

    forAll (pointsOnCellsToRefine, pointI)
    {
        if (pointsOnCellsToRefine[pointI])
        {
            // This point is an unrefinement candidate, set it and increment
            unrefinementCandidates[n++] = pointI;
        }
    }







//    // Loop through all split points and select candidates to unrefine
//    forAll (pField, pointI)
//    {
//        // Get point value
//        const scalar pointValue = pFieldIn[pointI];

//        if (pointValue < -0.5)
//        {
//            // Point value is outside of bounds, append point for potential
//            // unrefinement
//            unrefinementCandidates.append(pointI);
//        }

////        if
////        (
//////            pointValue < upperConstrainBound_
////            pointValue > lowerConstrainBound_
////        )
////        {
////            excludedUnrefCandidates.append(pointI);
////        }
////        else
////        {
////            // Point value is not within the intervall 
////            // where unrefinement is constrained
////            if (pointValue > upperBound_ || pointValue < lowerBound_)
//////            if (pointValue < (lowerBound_-1))
////            {
////                // Point value is outside of bounds, append point for potential
////                // unrefinement
////                unrefinementCandidates.append(pointI);
////            }
////        }
//    }

    // Print out some information
    Info<< "Selection algorithm " << type() << " selected "
        << returnReduce(unrefinementCandidates.size(), sumOp<label>())
        << " points as unrefinement candidates." 
        << endl;

    // Return the list in the Xfer container to prevent copying
    return unrefinementCandidates.xfer();
}


// ************************************************************************* //
