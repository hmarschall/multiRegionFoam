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
    Constantin Habes, TU Darmstadt.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fieldCurlRefinement.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "volPointInterpolation.H"
#include "pointVolInterpolation.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(fieldCurlRefinement, 0);
addToRunTimeSelectionTable
(
    refinementSelection,
    fieldCurlRefinement,
    dictionary
);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldCurlRefinement::fieldCurlRefinement
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    refinementSelection(mesh, dict),
    fieldName_(coeffDict().lookup("fieldName")),
    lowerBound_(readScalar(coeffDict().lookup("lowerBound"))),
    upperBound_(readScalar(coeffDict().lookup("upperBound"))),
    cellPointCellSmoothing_
    (
        coeffDict().lookupOrDefault<Switch>("cellPointCellSmoothing", false)
    )
{}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::fieldCurlRefinement::~fieldCurlRefinement()
{}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

Foam::Xfer<Foam::labelList>
Foam::fieldCurlRefinement::refinementCellCandidates() const
{
    // Create storage for collection of cells. Assume that one in five cells
    // will be refined to prevent excessive resizing.
    dynamicLabelList refinementCandidates(Foam::max(100, mesh().nCells()/5));

    // Get the vector field
    const volVectorField& vField = mesh().lookupObject<volVectorField>(fieldName_);

    // Create temporary for the curl magnitute field
    tmp<volScalarField> tempCurlMagFld(Foam::mag(fvc::curl(vField)));

    // Get const reference to curl magnitute field
    const volScalarField& curlMagFld = tempCurlMagFld();

    // Loop through internal field and collect cells to refine
    const scalarField& curlMagInFld = curlMagFld.internalField();

    forAll (curlMagInFld, cellI)
    {
        if (curlMagInFld[cellI] > lowerBound_ && curlMagInFld[cellI] < upperBound_)
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
Foam::fieldCurlRefinement::unrefinementPointCandidates() const
{
    // Create storage for collection of candidates. Assume that one in ten
    // mesh points will be unrefined to prevent excessive resizing
    dynamicLabelList unrefinementCandidates(Foam::max(100, mesh().nPoints()/10));

    // Get the vector field
    const volVectorField& vField = mesh().lookupObject<volVectorField>(fieldName_);

    // Create temporary for the curl magnitute field
    tmp<volScalarField> tempCurlMagFld(Foam::mag(fvc::curl(vField)));

    // Get const reference to curl magnitute field
    const volScalarField& curlMagFld = tempCurlMagFld();

    // Create volume to point interpolation object
    const volPointInterpolation& vpi = volPointInterpolation::New(mesh());

    // Interpolate the volume field from cell centres to points and get the
    // internal point field
    pointScalarField pField(vpi.interpolate(curlMagFld));
    const scalarField& pFieldIn = pField.internalField();

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
