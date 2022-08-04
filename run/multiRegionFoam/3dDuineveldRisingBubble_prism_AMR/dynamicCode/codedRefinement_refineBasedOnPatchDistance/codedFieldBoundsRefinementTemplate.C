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
    Mohammed Elwardi Fadeli (elwardifadeli@gmail.com)

\*---------------------------------------------------------------------------*/

#include "codedFieldBoundsRefinementTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "volPointInterpolation.H"
#include "pointVolInterpolation.H"

// Include user-supplied headers
#line 153 "::refineBasedOnPatchDistance::refinementSelection"
#include "volFields.H"
                    #include "patchWave.H"
                    #include "polyPatchID.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // A unique Function to make sure the correct library is loaded
    void codedRefinement_refineBasedOnPatchDistance_b30975726288a68b1a09207309e0b5f088bb938d(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}


defineTypeNameAndDebug(codedRefinement_refineBasedOnPatchDistance, 0);
addToRunTimeSelectionTable
(
    refinementSelection,
    codedRefinement_refineBasedOnPatchDistance,
    dictionary
);

// SHA1 as a static member variable
const char* const codedRefinement_refineBasedOnPatchDistance::SHA1sum =
    "b30975726288a68b1a09207309e0b5f088bb938d";

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::codedRefinement_refineBasedOnPatchDistance::codedRefinement_refineBasedOnPatchDistance
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    refinementSelection(mesh, dict),
    fieldName_(coeffDict().lookup("fieldName")),
    field_
    (
        IOobject
        (
            fieldName_,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    lowerBound_(readScalar(coeffDict().lookup("lowerBound"))),
    upperBound_(readScalar(coeffDict().lookup("upperBound"))),
    cellPointCellSmoothing_
    (
        coeffDict().lookupOrDefault<Switch>("cellPointCellSmoothing", false)
    )
{
    if (false)
    {
        Info<<"construct codedRefinement_refineBasedOnPatchDistance sha1: b30975726288a68b1a09207309e0b5f088bb938d from mesh\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::codedRefinement_refineBasedOnPatchDistance::~codedRefinement_refineBasedOnPatchDistance()
{
    if (false)
    {
        Info<<"destroy codedRefinement_refineBasedOnPatchDistance sha1: b30975726288a68b1a09207309e0b5f088bb938d\n";
    }
}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

Foam::Xfer<Foam::labelList>
Foam::codedRefinement_refineBasedOnPatchDistance::refinementCellCandidates() const
{
    // Do user calculations here
    {
        #line 160 "::refineBasedOnPatchDistance::refinementSelection"
Info << "*** dynamicCode ***" << endl;

                    // Read patch names from dictionary
                    const wordList patchNames(coeffDict().lookup("distancePatches"));

                    // Collect patchIDs into a hash set
                    labelHashSet patchIDs(patchNames.size());
                    forAll (patchNames, patchI)
                    {
                        // Get polyPatchID
                        const polyPatchID pID(patchNames[patchI], mesh().boundaryMesh());

                        if (pID.active())
                        {
                            patchIDs.insert(pID.index());
                        }
                        else
                        {
                            FatalIOErrorIn
                            (
                                "Xfer<labelList> minPatchDistanceRefinement::"
                                "refinenementPatchCandidates() const",
                                coeffDict()
                            )   << "Cannot find patch " << patchNames[patchI]
                                << " in the mesh." << nl
                                << "Available patches are: " << mesh().boundaryMesh().names()
                                << abort(FatalIOError);
                        }
                    }

                    patchWave waveDistance(mesh(), patchIDs, false);
                    const scalarField& patchDistance = waveDistance.distance();

                    scalar maxPatchDistance = gMax(patchDistance);
                    scalar minPatchDistance = gMin(patchDistance);

                    Info << "max patchDistance: " << maxPatchDistance << endl;
                    Info << "min patchDistance: " << minPatchDistance << endl;
                    
                    field_.internalField() = patchDistance;

                    field_.write();
                    mesh().write();
    }

    // Get the field
    const volScalarField& vField = field_;

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

        // bugfix: Avoid client interpolation at processor patches
        forAll(mesh().boundaryMesh(), pi)
        {
            const polyPatch& patch = mesh().boundaryMesh()[pi];
            if (isA<processorPolyPatch>(patch))
            {
                const labelList& faceCells = patch.faceCells();
                forAll(patch.meshPoints(), ni)
                {
                    const label& point = patch.meshPoints()[ni];
                    const labelList& faces = patch.pointFaces()[ni];
                    label acceptedFieldValue = vField[faceCells[faces[0]]];
                    forAll(faces, fi)
                    {
                        scalar value = vField[faceCells[faces[fi]]];
                        bool toBeConsidered = (value > lowerBound_ &&  value < upperBound_);
                        if (!toBeConsidered)
                        {
                            acceptedFieldValue = value;
                        }
                    }
                    pField[point] = acceptedFieldValue;
                }
            }
        }

        // Assign to the temporary object
        tvf = tInterpolatedVf;
    }

    // Get const reference to (possibly smoothed volume field)
    const volScalarField& vf = tvf();

    // Create storage for collection of cells. Assume that one in five cells
    // will be refined to prevent excessive resizing.
    dynamicLabelList refinementCandidates(Foam::max(100, mesh().nCells()/5));

    // Loop through internal field and collect cells to refine
    const scalarField& vfIn = vf.internalField();

    forAll (vfIn, cellI)
    {
        // Get current cell value
        const scalar& cellValue = vfIn[cellI];

        if (cellValue > lowerBound_ && cellValue < upperBound_)
        {
            // Cell value is within the bounds, append cell for potential
            // refinement
            refinementCandidates.append(cellI);
        }
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
Foam::codedRefinement_refineBasedOnPatchDistance::unrefinementPointCandidates() const
{
    // Do optional user calculations here
    {
        
    }

    // Get the field
    const volScalarField& vField = field_;

    // Create volume to point interpolation object
    const volPointInterpolation& vpi = volPointInterpolation::New(mesh());

    // Interpolate the volume field from cell centres to points and get the
    // internal point field
    pointScalarField pField(vpi.interpolate(vField));
    const scalarField& pFieldIn = pField.internalField();

    // bugfix: Avoid client interpolation at processor patches
    forAll(mesh().boundaryMesh(), pi)
    {
        const polyPatch& patch = mesh().boundaryMesh()[pi];
        if (isA<processorPolyPatch>(patch))
        {
            const labelList& faceCells = patch.faceCells();
            forAll(patch.meshPoints(), ni)
            {
                const label& point = patch.meshPoints()[ni];
                const labelList& faces = patch.pointFaces()[ni];
                label acceptedFieldValue = vField[faceCells[faces[0]]];
                forAll(faces, fi)
                {
                    scalar value = vField[faceCells[faces[fi]]];
                    bool toBeConsidered = (value > upperBound_ ||  value < lowerBound_);
                    if (toBeConsidered)
                    {
                        acceptedFieldValue = value;
                    }
                }
                pField[point] = acceptedFieldValue;            }
        }
    }

    // Create storage for collection of candidates. Assume that one in ten
    // mesh points will be unrefined to prevent excessive resizing
    dynamicLabelList unrefinementCandidates
    (
        Foam::max(100, mesh().nPoints()/10)
    );

    // Loop through all split points and select candidates to unrefine
    forAll (pField, pointI)
    {
        // Get point value
        const scalar pointValue = pFieldIn[pointI];

        if (pointValue > upperBound_ || pointValue < lowerBound_)
        {
            // Point value is outside of bounds, append point for potential
            // unrefinement
            unrefinementCandidates.append(pointI);
        }
    }

    // Print out some information
    Info<< "Selection algorithm " << type() << " selected "
        << returnReduce(unrefinementCandidates.size(), sumOp<label>())
        << " points as unrefinement candidates."
        << endl;

    // Return the list in the Xfer container to prevent copying
    return unrefinementCandidates.xfer();
}


// ************************************************************************* //

