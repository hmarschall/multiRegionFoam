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

#include "interfaceTrackingFvMesh/movingInterfacePatches.H"
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::interfaceTrackingFvMesh::correctPointNormals()
{
    forAll (movingInterfaceEntries_, surfI)
    {
        word patchName = movingInterfaceEntries_[surfI].keyword();

        movingInterfaces_[patchName]->correctPointNormals();
    }
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
    active_(true),
    timeIndex_(0),
    doMeshMotion_(false)
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

    if (movingInterfaceEntries_.empty())
    {
        active_ = false;
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

bool Foam::interfaceTrackingFvMesh::movePointsSurfacePatches()
{
    return true;
}

bool Foam::interfaceTrackingFvMesh::update()
{

    this->updatePatchMotionAndMoveMesh();

    if (debug)
    {
        const objectRegistry& db = this->thisDb();
        HashTable<const meshObjectBase*> tbl =
            db.lookupClass<meshObjectBase>();

        Info << nl <<"Mesh Objects registered to mesh :" << nl << endl;

        for
        (
            HashTable<const meshObjectBase*>::iterator iter =
                tbl.begin();
            iter != tbl.end();
            ++iter
        )
        {
            const meshObjectBase& obj = *(iter());

            Info << " type : " << obj.type() << nl
                 << endl;
        }
    }

    // Correct point normals
    //correctPointNormals();

	// Move correctedFvPatchField fvSubMeshes
	forAll(U().boundaryField(), patchI)
	{
		if
		(
			(
				U().boundaryField()[patchI].type()
			 == fixedGradientCorrectedFvPatchField<vector>::typeName
			)
			||
			(
				U().boundaryField()[patchI].type()
			 == fixedValueCorrectedFvPatchField<vector>::typeName
			)
			||
			(
				U().boundaryField()[patchI].type()
			 == zeroGradientCorrectedFvPatchField<vector>::typeName
			)
		)
		{
			correctedFvPatchField<vector>& aU =
				refCast<correctedFvPatchField<vector> >
				(
					const_cast<volVectorField&>(U()).boundaryField()[patchI]
				);

			aU.movePatchSubMesh();
		}
	}

	forAll(p().boundaryField(), patchI)
	{
		if
		(
			(
				p().boundaryField()[patchI].type()
			 == fixedGradientCorrectedFvPatchField<scalar>::typeName
			)
			||
			(
				p().boundaryField()[patchI].type()
			 == fixedValueCorrectedFvPatchField<scalar>::typeName
			)
			||
			(
				p().boundaryField()[patchI].type()
			 == zeroGradientCorrectedFvPatchField<scalar>::typeName
			)
		)
		{
			correctedFvPatchField<scalar>& aP =
				refCast<correctedFvPatchField<scalar> >
				(
					const_cast<volScalarField&>(p()).boundaryField()[patchI]
				);

			aP.movePatchSubMesh();
		}
	}

    // Mesh motion only - return true since flux re-calculation needed
    return true;
}

void Foam::interfaceTrackingFvMesh::updatePatchMotionAndMoveMesh()
{
    this->updatePatchMotion();

    Info << nl << "Updating " << this->name() << nl
         << "Perform " << (doMeshMotion_ ? "" : "no ")
         << "whole mesh motion " << nl << endl;

    if (doMeshMotion_)
    {
            motionPtr_->solve();

            fvMesh::movePoints(motionPtr_->curPoints());

            forAll (movingInterfaceEntries_, surfI)
            {
                word patchName = movingInterfaceEntries_[surfI].keyword();

                if
                (
                    !movingInterfaces_[patchName]->isFsiInterface()
                )
                {
                    movingInterfaces_[patchName]->updateDisplacementDirections();
                }
            }

            // Reset doMeshMotion
            doMeshMotion_ = false;
    }
}

void Foam::interfaceTrackingFvMesh::updatePatchMotion()
{
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

            // switch true as soon as one of the moving patches
            // DOES trigger mesh motion
            doMeshMotion_ =
            (
                doMeshMotion_
             || movingInterfaces_[patchName]->moveMesh()
            );

            if
            (
                movingInterfaces_[patchName]->isInterface()
                && !movingInterfaces_[patchName]->isFsiInterface()
            )
            {
                refCast<const interfaceTrackingFvMesh>
                (
                    movingInterfaces_[patchName]->nbrMesh()
                ).doMeshMotion(doMeshMotion_);
            }

            if (movingInterfaces_[patchName]->moveMesh())
            {
                motionUPatch ==
                    surfacePointDisplacement/mesh().time().deltaT().value();
            }
            else
            {
                vectorField zeroDisplacement
                    (
                        movingInterfaces_[patchName]->
                            aMesh().patch().localPoints().size(),
                        vector::zero
                    );

                motionUPatch ==
                    zeroDisplacement/mesh().time().deltaT().value();
            }

            // Set shadow patch motion (interface)
            if
            (
                movingInterfaces_[patchName]->isInterface()
                && !movingInterfaces_[patchName]->isFsiInterface()
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

                if(movingInterfaces_[patchName]->moveMesh())
                {
                    motionNbrUPatch ==
                        movingInterfaces_[patchName]->shadowPointDisplacement
                        (
                            surfacePointDisplacement
                        )/mesh().time().deltaT().value();

                }
                else
                {
                    vectorField zeroDisplacement
                        (
                            movingInterfaces_[patchName]->
                                aMesh().patch().localPoints().size(),
                            vector::zero
                        );

                    motionNbrUPatch ==
                        movingInterfaces_[patchName]->shadowPointDisplacement
                        (
                            zeroDisplacement
                        )/mesh().time().deltaT().value();
                }
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

            //correctPointNormals();

            // switch true as soon as one of the moving patches
            // DOES trigger mesh motion
            doMeshMotion_ =
            (
                doMeshMotion_
             || movingInterfaces_[patchName]->moveMesh()
            );

            if
            (
                movingInterfaces_[patchName]->isInterface()
                && !movingInterfaces_[patchName]->isFsiInterface()
            )
            {
                refCast<const interfaceTrackingFvMesh>
                (
                    movingInterfaces_[patchName]->nbrMesh()
                ).doMeshMotion(doMeshMotion_);
            }

            if (movingInterfaces_[patchName]->moveMesh())
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
                tetPolyPatchInterpolation tppiPatch
                (
                    refCast<const faceTetPolyPatch>
                    (
                        motionUPatch.patch()
                    )
                );

                vectorField zeroDisplacement
                    (
                        movingInterfaces_[patchName]->
                            aMesh().patch().localPoints().size(),
                        vector::zero
                    );

                motionUPatch ==
                    tppiPatch.pointToPointInterpolate
                    (
                        zeroDisplacement/mesh().time().deltaT().value()
                    );
            }

            // Set shadow patch motion (interface)
            if
            (
                movingInterfaces_[patchName]->isInterface()
                && !movingInterfaces_[patchName]->isFsiInterface()
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


                if(movingInterfaces_[patchName]->moveMesh())
                {
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
                else
                {
                    tetPolyPatchInterpolation tppiPatch
                    (
                        refCast<const faceTetPolyPatch>
                        (
                            motionNbrUPatch.patch()
                        )
                    );

                    vectorField zeroDisplacement
                    (
                        movingInterfaces_[patchName]->
                            aMesh().patch().localPoints().size(),
                        vector::zero
                    );

                    motionNbrUPatch ==
                        tppiPatch.pointToPointInterpolate
                        (
                            movingInterfaces_[patchName]->shadowPointDisplacement
                            (
                                zeroDisplacement
                            )/mesh().time().deltaT().value()
                        );

                }
            }
        }

        Info<< "Patch motion trigger "
            << patchName << ": "
            << movingInterfaces_[patchName]->movePatch() << endl;

        if (movingInterfaces_[patchName]->movePatch())
        {
            movingInterfaces_[patchName]->enforcePatchMotion();

            if
            (
                movingInterfaces_[patchName]->isInterface()
                && !movingInterfaces_[patchName]->isFsiInterface()
            )
            {
                movingInterfaces_[patchName]->enforceShadowPatchMotion();
            }

            // Correct point normals
            //correctPointNormals();
        }
    }
}


// ************************************************************************* //
