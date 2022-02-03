/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Zeljko Tukovic, FSB Zagreb.
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "interfaceTrackingFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "volFields.H"

#include "velocityLaplacianFvMotionSolver.H"
#include "twoDPointCorrector.H"
#include "demandDrivenData.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interfaceTrackingFvMesh, 0);
    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        interfaceTrackingFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceTrackingFvMesh::interfaceTrackingFvMesh
(
    const IOobject& io
)
:
    topoChangerFvMesh(io),
    motionPtr_(motionSolver::New(*this)),
    motionDict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                io.time().constant(),
                (io.name() == topoChangerFvMesh::defaultRegion ? "" : io.name() ),
                io.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).subDict(typeName + "Coeffs")
    ),
    movingInterfaces_(),
    movingInterfaceEntries_
    (
        motionDict_.lookup("movingSurfacePatches")
    )
{
//    movingInterfaces_.setSize(movingInterfaceEntries_.size());

    forAll (movingInterfaceEntries_, surfI)
    {
        movingInterfaces_.set
        (
//            surfI,
            movingInterfaceEntries_[surfI].keyword(),
            new movingInterfacePatches
            (
                movingInterfaceEntries_[surfI].keyword(),
                *this,
                movingInterfaceEntries_[surfI].dict()
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceTrackingFvMesh::~interfaceTrackingFvMesh()
{
    this->clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::volVectorField& Foam::interfaceTrackingFvMesh::U() const
{
    return mesh().objectRegistry::lookupObject<volVectorField>("U");
}


const Foam::volScalarField& Foam::interfaceTrackingFvMesh::p() const
{
    return mesh().objectRegistry::lookupObject<volScalarField>("p");
}


const Foam::surfaceScalarField& Foam::interfaceTrackingFvMesh::phi() const
{
    return mesh().objectRegistry::lookupObject<surfaceScalarField>("phi");
}


bool Foam::interfaceTrackingFvMesh::update()
{
    Info << nl << "*** Updating " << this->name() << endl;

//        fvMotionSolver& mSolver =
//            dynamic_cast<fvMotionSolver&>
//            (
//                motionPtr_()
//            );

//        pointVectorField& motionU = mSolver.pointMotionU();

    forAll (movingInterfaceEntries_, surfI)
    {
        // Set surface patch motion
        word patchName = movingInterfaceEntries_[surfI].keyword();

        // Update
//        movingInterfaces_[patchName]->updateInterpolatorAndGlobalPatches();

        pointVectorField& motionU =
            const_cast<pointVectorField&>
            (
                mesh().objectRegistry::
                lookupObject<pointVectorField>
                (
                    "pointMotionU"
                )
            );

        fixedValuePointPatchVectorField& motionUPatch =
            refCast<fixedValuePointPatchVectorField>
            (
                motionU.boundaryField()
                [movingInterfaces_[patchName]->patchID()]
            );

        //motionUPatch ==
        //    totalDisplacement()/mesh().time().deltaT().value();

//        motionUPatch ==
//            movingInterfaces_[patchName]->surfacePointDisplacement()
//            /mesh().time().deltaT().value();

        vectorField surfacePointDisplacement
            (movingInterfaces_[patchName]->surfacePointDisplacement());

        motionUPatch ==
            surfacePointDisplacement/mesh().time().deltaT().value();

        // Set shadow patch motion (interface)
        if (movingInterfaces_[patchName]->isInterface())
        {
            Info << nl << "*** Moving points of interface shadow" << endl;

            pointVectorField& motionNbrU =
                const_cast<pointVectorField&>
                (
                    movingInterfaces_[patchName]->nbrMesh().objectRegistry::
                    lookupObject<pointVectorField>
                    (
                        "pointMotionU"
                    )
                );

            fixedValuePointPatchVectorField& motionNbrUPatch =
                refCast<fixedValuePointPatchVectorField>
                (
                    motionNbrU.boundaryField()
                    [movingInterfaces_[patchName]->nbrPatch().index()]
                );

            motionNbrUPatch ==
                movingInterfaces_[patchName]->shadowPointDisplacement
                (
                    surfacePointDisplacement
                )/mesh().time().deltaT().value();
        }
    }

//    twoDPointCorrector twoDPointCorr(mesh());
//    twoDPointCorr.correctPoints(newMeshPoints);

//    fvMesh::movePoints(newMeshPoints);
//    setOldPoints(newMeshPoints);
//    movePoints(newMeshPoints);

//    aMesh().movePoints();

    motionPtr_->solve();
//-    mSolver.solve();

//    dynamicMotionSolverFvMesh::update();

    fvMesh::movePoints(motionPtr_->curPoints());
//-    fvMesh::movePoints(mSolver.curPoints());

//    fvMesh::movePoints(mSolver.newPoints());

    return true;
}


// ************************************************************************* //
