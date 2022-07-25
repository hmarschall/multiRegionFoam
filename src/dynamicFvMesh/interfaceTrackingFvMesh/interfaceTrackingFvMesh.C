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
#include "velocityLaplacianFvMotionSolver.H"
#include "laplaceTetMotionSolver.H"

#include "tetPolyPatchFields.H"
#include "fixedValueTetPolyPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "tetFemMatrices.H"
#include "tetPointFields.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"

#include "slipFvPatchFields.H"
#include "wedgeFvsPatchFields.H"
#include "symmetryFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "zeroGradientCorrectedFvPatchFields.H"
#include "fixedGradientCorrectedFvPatchFields.H"
#include "fixedValueCorrectedFvPatchFields.H"

#include "volFields.H"
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
    fvMotionSolver_(false),
    feMotionSolver_(false),
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
    ),
    active_(true)
{
    // Set movingInterfacePatches objects
    forAll (movingInterfaceEntries_, surfI)
    {
        movingInterfaces_.set
        (
            movingInterfaceEntries_[surfI].keyword(),
            new movingInterfacePatches
            (
                movingInterfaceEntries_[surfI].keyword(),
                *this,
                movingInterfaceEntries_[surfI].dict()
            )
        );
    }

    Info << "Motion solver type : " << motionPtr_->type() << endl;

    // Check mesh motion solver type
    if
    (
        motionPtr_->type()
     == laplaceTetMotionSolver::typeName
    )
    {
        feMotionSolver_ =
            mesh().objectRegistry::foundObject<tetPointVectorField>
            (
                "motionU"
            );
    }
    else if
    (
        motionPtr_->type()
     == velocityLaplacianFvMotionSolver::typeName
    )
    {
        fvMotionSolver_ =
            mesh().objectRegistry::foundObject<pointVectorField>
            (
                "pointMotionU"
            );
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported mesh motion solver : "
            << motionPtr_->type() << nl
            << abort(FatalError);
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
//    if (!active_)
//    {
//        return true;
//    }

//    bool control = false;

    forAll (movingInterfaceEntries_, surfI)
    {
        // Set surface patch motion
        word patchName = movingInterfaceEntries_[surfI].keyword();

        if (fvMotionSolver_)
        {
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

            vectorField surfacePointDisplacement
                (movingInterfaces_[patchName]->surfacePointDisplacement());

            surfacePointDisplacement += 
                movingInterfaces_[patchName]->smoothSurfaceMesh();

            motionUPatch ==
                surfacePointDisplacement/mesh().time().deltaT().value();

            // TODO: needs thoughts on logics for multi-region arrangement
//            dynamicFvMesh& aleNbgMesh =
//                const_cast<dynamicFvMesh&>
//                (
//                    movingInterfaces_[patchName]->nbrMesh()
//                );

//            interfaceTrackingFvMesh& aleNbrMesh =
//                refCast<interfaceTrackingFvMesh>
//                (
//                    aleNbgMesh
//                );

//            if (movingInterfaces_[patchName]->moving())
//            {
//                motionUPatch ==
//                    surfacePointDisplacement/mesh().time().deltaT().value();

//                if (!aleNbrMesh.active())
//                {
//                    control = true;

//                    aleNbrMesh.active() = true;
//                }
//            }
//            else
//            {
//                if (!control)
//                {
//                    aleNbrMesh.active() = false;
//                }

//                return true;
//            }

            // Set shadow patch motion (interface)
            if 
            (
                movingInterfaces_[patchName]->isInterface()
             && movingInterfaces_[patchName]->moving()
            )
            {
                pointVectorField& motionNbrU =
                    const_cast<pointVectorField&>
                    (
                        movingInterfaces_[patchName]->nbrMesh()
                        .objectRegistry::lookupObject<pointVectorField>
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

        if (feMotionSolver_)
        {
            tetPointVectorField& motionU =
                const_cast<tetPointVectorField&>
                (
                    mesh().objectRegistry::
                    lookupObject<tetPointVectorField>
                    (
                        "motionU"
                    )
                );

            fixedValueTetPolyPatchVectorField& motionUPatch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()
                    [movingInterfaces_[patchName]->patchID()]
                );

            vectorField surfacePointDisplacement
                (movingInterfaces_[patchName]->surfacePointDisplacement());

            surfacePointDisplacement += 
                movingInterfaces_[patchName]->smoothSurfaceMesh();

            if (movingInterfaces_[patchName]->moving())
            {
                tetPolyPatchInterpolation tppiPatch
                (
                    refCast<const faceTetPolyPatch>
                    (
                        motionUPatch.patch()
                    )
                );

                motionUPatch ==
                    tppiPatch.pointToPointInterpolate
                    (
                        surfacePointDisplacement/mesh().time().deltaT().value()
                    );
            }
            else
            {
                return true;
            }

            // Set shadow patch motion (interface)
            if 
            (
                movingInterfaces_[patchName]->isInterface()
             && movingInterfaces_[patchName]->moving()
            )
            {
                tetPointVectorField& motionNbrU =
                    const_cast<tetPointVectorField&>
                    (
                        movingInterfaces_[patchName]->nbrMesh()
                        .objectRegistry::lookupObject<tetPointVectorField>
                        (
                            "motionU"
                        )
                    );

                fixedValueTetPolyPatchVectorField& motionNbrUPatch =
                    refCast<fixedValueTetPolyPatchVectorField>
                    (
                        motionNbrU.boundaryField()
                        [movingInterfaces_[patchName]->nbrPatch().index()]
                    );

                tetPolyPatchInterpolation tppiPatch
                (
                    refCast<const faceTetPolyPatch>
                    (
                        motionNbrUPatch.patch()
                    )
                );

                motionNbrUPatch ==
                    tppiPatch.pointToPointInterpolate
                    (
                        movingInterfaces_[patchName]->shadowPointDisplacement
                        (
                            surfacePointDisplacement
                        )/mesh().time().deltaT().value()
                    );
            }
        }
    }

    Info << nl << "Updating " << this->name() << nl << endl;

    motionPtr_->solve();

    fvMesh::movePoints(motionPtr_->curPoints());

    // Correct point normals and finite area mesh
    forAll (movingInterfaceEntries_, surfI)
    {
        word patchName = movingInterfaceEntries_[surfI].keyword();

        movingInterfaces_[patchName]->correctPointNormals();

        movingInterfaces_[patchName]->aMesh().movePoints();
    }

    return true;
}


// ************************************************************************* //
