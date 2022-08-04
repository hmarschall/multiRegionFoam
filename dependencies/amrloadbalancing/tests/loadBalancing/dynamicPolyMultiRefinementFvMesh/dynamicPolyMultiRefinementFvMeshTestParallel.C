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

Description
    Parallel Tests for dynamicPolyMultiRefinementFvMesh class

Author
    Mohammed Elwardi Fadeli (elwardifadeli@gmail.com)

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "catch.H"
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "polyTopoChanger.H"
#include "SortableList.H"
#include "IStringStream.H"
#include "dynamicPolyMultiRefinementFvMesh.H"
#include "ReadFields.H"
#include "IOobjectList.H"
#include "labelList.H"
#include "boxToCell.H"
#include "cellSet.H"
#include "memberStealer.H"
#include "testMacros.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace multiCriteriaParallel
{
    // Requirement 0-0: Time for parallel runs
    #include "createTestTime.H"

    // Requirement 0-1: Pointer for mesh
    Foam::autoPtr<Foam::dynamicPolyMultiRefinementFvMesh> meshPtr;

} // End of namespace multiCriteriaParallel

using namespace Foam;
using namespace multiCriteriaParallel;

prepareTimePaths();

SCENARIO
(
    "Multi-criteria adaptive poly refinement in parallel",
    "[Integration][Parallel][!throws]"
)
{
    // Switch to processor case
    word newCase = "processor"+Foam::name(Pstream::myProcNo());
    modifyTimePaths(runTime, true, newCase);

    // Test with 8 different refinement levels
    auto maxRefLevel1 = GENERATE( range(1,3) );
    auto maxRefLevel2 = GENERATE( range(1,3) );
    
    // Stuff to capture for meaningfull error messages
    CAPTURE
    (
        maxRefLevel1,
        maxRefLevel2,
        Pstream::myProcNo(),
        runTime.caseName()
    );

    GIVEN
    (
        "A time object, and a dynamicMeshDict on processor "
        + Foam::name(Pstream::myProcNo())
    )
    {
        // Adds coeffs subdict with 2 refinement selections
        // 1st: static box (alpha0)
        // 2nd: moving interface (alpha1)
        IStringStream is
        (
            "dynamicRefineFvMeshCoeffs{                                    "
            "maxRefinementLevel   1;                                       "
            "refineInterval   1;                                           "
            "unrefineInterval 1;                                           "
            "separateUpdates false;                                        "
            "active yes;                                                   "
            "maxCells             2000000;                                 "
            "nRefinementBufferLayers        1;                             "
            "nUnrefinementBufferLayers      3;                             "
            "edgeBasedConsistency           yes;                           "
            "refinements                                                   "
            "(                                                             "
            "    basedOnAlpha0                                             "
            "    {                                                         "
            "        maxRefinementLevel   "+Foam::name(maxRefLevel1)+";    "
            "        refinementSelection                                   "
            "        {                                                     "
            "            type       fieldBoundsRefinement;                 "
            "            fieldName  alpha0;                                "
            "            lowerBound 0.499;                                 "
            "            upperBound 0.501;                                 "
            "            cellPointCellSmoothing on;                        "
            "        }                                                     "
            "    }                                                         "
            "    basedOnAlpha1                                             "
            "    {                                                         "
            "        maxRefinementLevel   "+Foam::name(maxRefLevel2)+";    "
            "        refinementSelection                                   "
            "        {                                                     "
            "            type        interfaceRefinement;                  "
            "            fieldNames  ( alpha1 );                           "
            "            innerRefLayers  1;                                "
            "            outerRefLayers  1;                                "
            "            cellPointCellSmoothing on;                        "
            "        }                                                     "
            "    }                                                         "
            ");                                                            "
            "}                                                             "
        );

        IOdictionary dynamicMeshDict
        (
            IOobject
            (
                "dynamicMeshDict",
                runTime.constant(),
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            is
        );
        dynamicMeshDict.regIOobject::write();
        
        WHEN
        (
            "Two refinement criterions are applied, "
            " with maxRefinementLevels (" << maxRefLevel1
            << "," << maxRefLevel2 << ") "
            "each based on its own field"
        )
        {
            // Create the mesh with default name
            resetMeshPointer
            (
                runTime,
                meshPtr,
                dynamicPolyMultiRefinementFvMesh,
                dynamicFvMesh::defaultRegion
            );
            auto& mesh = meshPtr();

            // Create fields to control refinement
            createField(alpha0, volScalarField, runTime, mesh, dimless, 0.0);
            createField(alpha1, volScalarField, runTime, mesh, dimless, 0.0);

            // Initialize field values
            // (minX maxX minY maxY)
            std::vector<scalar> boxBounds = {0.2, 0.4, 0.25, 0.35};
            std::vector<scalar> interfaceBounds = {0.1, 0.12};
            forAll(mesh.C(), ci)
            {
                const auto& center = mesh.C()[ci];
                // Set alpha0 to 0.5 inside a box
                if 
                (
                    center.x() >= boxBounds[0] && center.x() <= boxBounds[1]
                    && center.y() >= boxBounds[2] && center.y() <= boxBounds[3]
                )
                alpha0[ci] = 0.5;
                // Set alpha1 to 1 at some fixed x-coordinate
                if 
                (
                    center.x() >= interfaceBounds[0]
                    && center.x() <= interfaceBounds[1]
                )
                alpha1[ci] = 1.0;
            }

            // boxToCell configuration
            IStringStream boxIs("(0.2 0.25 -10) (0.4 0.35 10)");
            IStringStream interfaceIs("(0.1 -10 -10) (0.12 10 10)");

            // Initial cell count in refined box
            boxToCell box(mesh, boxIs);
            cellSet boxCells(mesh, "boxCells", IOobject::NO_READ);
            box.applyToSet(topoSetSource::setAction::NEW, boxCells);

            // Initial cell count in initial refined interface region
            boxToCell interface(mesh, interfaceIs);
            cellSet interfaceCells(mesh, "interfaceCells", IOobject::NO_READ);
            interface.applyToSet(topoSetSource::setAction::NEW, interfaceCells);

            bool boxDidntUnrefine = true;
            bool interfaceRefined = true;
            bool interfaceUnrefined = true;

            auto previousBoxCells = boxCells.size();

            // Simulate a time loop, this is typically done later in a THEN block;
            // but we want to test more than one aspect
            runTime.setTime(0.0, 0);
            while (runTime.run()) {
                runTime++;
                // Move alpha1 wave by 0.2 in X direction
                forAll(mesh.C(), ci)
                {
                    const auto& center = mesh.C()[ci];
                    label i = runTime.timeIndex() > 0
                        ? runTime.timeIndex() - 1 : 0;
                    if
                    (
                        center.x() >= interfaceBounds[0]+0.2*i
                        && center.x() <= interfaceBounds[1]+0.2*i)
                    {
                        alpha1[ci] = 1.0;
                    } else {
                        alpha1[ci] = 0.0;
                    }
                }

                // Update mesh
                mesh.update();
                
                // Write to disk
                runTime.writeNow();

                // Cell count in refined box
                IStringStream newBoxIs("(0.2 0.25 -10) (0.4 0.35 10)");
                boxToCell newBox(mesh, newBoxIs);
                cellSet newBoxCells(mesh, "newBoxCells", IOobject::NO_READ);
                newBox.applyToSet(topoSetSource::setAction::NEW, newBoxCells);

                // After interface passes the box, see if cell count ever decreases
                if(runTime.timeIndex() > 2)
                {
                    if(newBoxCells.size() < previousBoxCells)
                        boxDidntUnrefine = false;
                    previousBoxCells = newBoxCells.size();
                }

                // Cell count in original interface region
                IStringStream newInterfaceIs("(0.1 -10 -10) (0.12 10 10)");
                boxToCell newInterface(mesh, newInterfaceIs);
                cellSet newInterfaceCells(mesh, "interfaceCells", IOobject::NO_READ);
                newInterface.applyToSet(topoSetSource::setAction::NEW, newInterfaceCells);

                if (runTime.timeIndex() == 1)
                {
                    // If timeStep 1 just ended, see if interface cells were refined
                    if(newInterfaceCells.size() <= interfaceCells.size() && newInterfaceCells.size())
                        interfaceRefined = false;
                } else if (runTime.timeIndex() >= 2)
                {
                    // if timeStep 2 or later, see if interface cells were unrefined
                    if(newInterfaceCells.size() != interfaceCells.size() && newInterfaceCells.size())
                        interfaceUnrefined = false;
                }
            }

            THEN("Interface refines, then unrefine previous interface region")
            {
                REQUIRE(interfaceRefined);
                REQUIRE(interfaceUnrefined);

                AND_THEN("Second refinement doesn't unrefine previouly refined cells")
                {
                    REQUIRE(boxDidntUnrefine);
                }
            }

            // Free mesh memory
            meshPtr->clear();
        }
    }
    Pstream::waitRequests();
}

#include "undefineTestMacros.H"

// ************************************************************************* //
