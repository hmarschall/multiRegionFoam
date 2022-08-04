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
    Daniel Deising, TU Darmstadt. All rights reserved.
    Constantin Habes, TU Darmstadt. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "interfaceRefinement.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "volPointInterpolation.H"
#include "pointVolInterpolation.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(interfaceRefinement, 0);
addToRunTimeSelectionTable
(
    refinementSelection,
    interfaceRefinement,
    dictionary
);

}

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

void Foam::interfaceRefinement::interfaceCandidates() const
{
    forAll (isInterface_, i)
    {
        isInterface_[i] *= 0.0;
    }

    forAll(fldNames_, n)
    {
        // Get the field
        const volScalarField& vField = 
            mesh().lookupObject<volScalarField>(fldNames_[n]);

        // Create temporary for the magnitute of the surface gradient field
        surfaceScalarField deltaAlpha = 
            mag(fvc::snGrad(vField)/mesh().deltaCoeffs());
            // /mag(gMax(vField)-gMin(vField));

        // List of owner und neighbour cells
        const unallocLabelList& owner = mesh().owner();
        const unallocLabelList& neighbour = mesh().neighbour();

        forAll(deltaAlpha, faceI)
        {
            // currently a fixed prescribed value; 
            // should be read from díctionary
            if (deltaAlpha[faceI] > tolerance_)
            {
                label own = owner[faceI];
                label nei = neighbour[faceI];

                // set isInterface field to one
                isInterface_[n][own] = 1.0;
                isInterface_[n][nei] = 1.0;
            }
        }

        // assumed 0.5*(vFieldmin+vFieldmax) defines the interface
        dimensionedScalar vFieldInterfaceValue
        (
            0.5*(gMax(vField)+gMin(vField))
        );

        // Face-wise interplation for 2D meshes since point interpolation
        // does not work on empty patches
        if (mesh().nGeometricD() == 2)
        {
            //-DD: implementation based on face interpolation
            //     which results in slower transport in diagonal direction
            // add inner refinement layers
            for(label i=0; i < innerRefLayers_; i++)
            {
                isInterface_[n] += neg
                (
                    -fvc::average(fvc::interpolate(isInterface_[n]))
                    *pos(vField - vFieldInterfaceValue)
                );

                isInterface_[n] = neg(-isInterface_[n]);
            }
            
            // add outer refinement layers
            for(label i=0; i < outerRefLayers_; i++)
            {
                isInterface_[n] += neg
                (
                    -fvc::average(fvc::interpolate(isInterface_[n]))
                    *pos(vFieldInterfaceValue - vField)
                );

                isInterface_[n] = neg(-isInterface_[n]);
            }
        }
        else 
        {
            //-DD: version using volPointInterpolation 
            // (direction independent buffer layer)
            const volPointInterpolation& pInterp = 
                volPointInterpolation::New(mesh());

            // add inner refinement layers
            for(label i=0; i < innerRefLayers_; i++)
            {
                volScalarField markInner
                (
                    isInterface_[n]*pos(vField - vFieldInterfaceValue)
                );

                pointScalarField markLayerP(pInterp.interpolate(markInner));

                forAll(mesh().C(), cellI)
                {
                    scalar sum = 0.;
                    label nPoints = 0;

                    forAll(mesh().cellPoints()[cellI], pointI)
                    {
                        sum += markLayerP[mesh().cellPoints()[cellI][pointI]];
                        nPoints++;
                    }
                    if (nPoints > 0)
                    {
                        sum /= nPoints;
                    }
                    isInterface_[n][cellI] += sum;
                }
            }

            isInterface_[n] = pos(isInterface_[n] - SMALL);

            // add outer refinement layers
            for(label i=0; i < outerRefLayers_; i++)
            {
                volScalarField markOuter
                (
                    isInterface_[n]*pos(vFieldInterfaceValue - vField)
                );

                pointScalarField markLayerP(pInterp.interpolate(markOuter));

                forAll(mesh().C(), cellI)
                {
                    scalar sum = 0.;
                    label nPoints = 0;

                    forAll(mesh().cellPoints()[cellI], pointI)
                    {
                        sum += markLayerP[mesh().cellPoints()[cellI][pointI]];
                        nPoints++;
                    }
                    if (nPoints > 0)
                    {
                        sum /= nPoints;
                    }
                    isInterface_[n][cellI] += sum;
                }
            }

            isInterface_[n] = pos(isInterface_[n] - SMALL);
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceRefinement::interfaceRefinement
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    refinementSelection(mesh, dict),
    fldNames_(coeffDict().lookup("fieldNames")),
    mesh_(mesh),
    innerRefLayers_(readScalar(coeffDict().lookup("innerRefLayers"))),
    outerRefLayers_(readScalar(coeffDict().lookup("outerRefLayers"))),
    cellPointCellSmoothing_
    (
        coeffDict().lookupOrDefault<Switch>("cellPointCellSmoothing", false)
    ),
    tolerance_(coeffDict().lookupOrDefault<scalar>("tolerance", 0.1)),
    isInterface_(fldNames_.size())
{
    forAll(fldNames_, i)
    {
        isInterface_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "isInterface" + Foam::name(i),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("isInterface", dimless, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::interfaceRefinement::~interfaceRefinement()
{}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

Foam::Xfer<Foam::labelList>
Foam::interfaceRefinement::refinementCellCandidates() const
{
    dynamicLabelList refinementCandidates
    (
        Foam::max(100, mesh().nCells()/5)
    );

    interfaceCandidates();

    forAll(isInterface_, i)
    {
        forAll(isInterface_[i], cellI)
        {
            if (isInterface_[i][cellI] > 0.5)
            {
                // Cell is around an interface
                // append to list of refinement candidates
                refinementCandidates.append(cellI);
            }
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
Foam::interfaceRefinement::unrefinementPointCandidates() const
{
    // Using volPointInterpolation of marker field
//    // Create storage for collection of cells. Assume that one in 10 cells
//    // will be unrefined to prevent excessive resizing.
//    dynamicLabelList unrefinementCandidates
//    (
//        Foam::max(100, mesh().nPoints()/10)
//    );

//    // Create volume to point interpolation object
//    const volPointInterpolation& vpi = volPointInterpolation::New(mesh());

//    // Interpolate the volume field from cell centres to points and
//    // get the internal point field
//    forAll (isInterface_, i)
//    {
//        pointScalarField pField(vpi.interpolate(isInterface_[i]));
//        const scalarField& pFieldIn = pField.internalField();

//        // Loop through all split points and select candidates to unrefine
//        forAll(pField, pointI)
//        {
//            if (pFieldIn[pointI] <= 0.5)
//            {
//                // Point is not around an interface
//                // append to list of unrefinement candidates
//                unrefinementCandidates.append(pointI);
//            }
//        }
//    }

    // Using cellPoints of marked cells
    boolList pointsOnInterfaceCells(mesh().nPoints(), false);

    // Get initial mesh data
    const labelListList& meshCellPoints = mesh().cellPoints();

    // Loop through all cells and count number of points of interface cells
    label nInterfacePoints = 0;

    forAll (isInterface_, i)
    {
        forAll (isInterface_[i], cellI)
        {
            const scalar& cInt = isInterface_[i].internalField()[cellI];

            if
            (
                (cInt > 0.5)
            )
            {
                // mark all cell points
                const labelList& curCellPoints = meshCellPoints[cellI];

                forAll (curCellPoints, n)
                {
                    // Get marker for this point
                    bool& ptOnRefCell = pointsOnInterfaceCells[curCellPoints[n]];

                    if (!ptOnRefCell)
                    {
                        // This points has not been marked yet, mark it and
                        // increment the counter for protected points
                        ptOnRefCell = true;
                        ++nInterfacePoints;
                    }
                }
            }
        }
    }

    // Create the list for unrefinement candidates
    labelList unrefinementCandidates(mesh().nPoints() - nInterfacePoints);
    label nUnrefPoints = 0;

    forAll (pointsOnInterfaceCells, pointI)
    {
        if (!pointsOnInterfaceCells[pointI])
        {
            // This point is an unrefinement candidate, set it and increment
            unrefinementCandidates[nUnrefPoints++] = pointI;
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
