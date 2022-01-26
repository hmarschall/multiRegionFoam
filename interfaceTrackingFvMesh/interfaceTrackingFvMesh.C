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
#include "wedgeFaPatch.H"
#include "wedgeFaPatchFields.H"
#include "slipFaPatchFields.H"
#include "fixedValueFaPatchFields.H"
#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"
#include "wallFvPatch.H"
#include "polyPatchID.H"
#include "fvcMeshPhi.H"
#include "velocityLaplacianFvMotionSolver.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"
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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::interfaceTrackingFvMesh::initializeData()
{
    // Set motion-based data
    movingSurfacePatches_ = wordList
    (
        motionDict_.lookup("movingSurfacePatches")
    );

    fixedSurfacePatches_ = wordList
    (
        motionDict_.lookup("fixedSurfacePatches")
    );

    pointNormalsCorrectionPatches_ = wordList
    (
        motionDict_.lookup("pointNormalsCorrectionPatches")
    );

    normalMotionDir_ = Switch
    (
        motionDict_.lookup("normalMotionDir")
    );

    smoothing_ = Switch
    (
        motionDict_.lookupOrDefault<Switch>("smoothing", false)
    );

    // Set surface patch index
    {
        const word surfacePatchName("freeSurface");

        polyPatchID patch(surfacePatchName, this->boundaryMesh());

        if (!patch.active())
        {
            FatalErrorInFunction
                << "Patch name " << surfacePatchName << " not found."
                << abort(FatalError);
        }

        surfacePatchID_ = patch.index();
    }

    // Set point normal correction for finite area mesh
    {
        boolList& correction = aMesh().correctPatchPointNormals();

        forAll(pointNormalsCorrectionPatches_, patchI)
        {
            word patchName = pointNormalsCorrectionPatches_[patchI];

            label patchID = aMesh().boundary().findPatchID(patchName);

            if(patchID == -1)
            {
                FatalErrorIn
                (
                    "surfaceTracking::surfaceTracking(...)"
                )   << "Patch name for point normals correction does not exist"
                    << abort(FatalError);
            }

            correction[patchID] = true;
        }
    }

    // Read motion direction
    if (!normalMotionDir_)
    {
        motionDir_ = vector(motionDict_.lookup("motionDir"));
        motionDir_ /= mag(motionDir_) + SMALL;
    }
}


void Foam::interfaceTrackingFvMesh::makeControlPoints()
{
    if (debug)
    {
        Info<< "interfaceTrackingFvMesh::makeControlPoints() : "
            << "making control points"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!controlPointsPtr_.empty())
    {
        FatalErrorIn("interfaceTrackingFvMesh::makeInterpolators()")
            << "patch to patch interpolators already exists"
            << abort(FatalError);
    }

    IOobject controlPointsHeader
    (
        "controlPoints",
        mesh().time().timeName(),
        mesh(),
        IOobject::MUST_READ
    );

    if (controlPointsHeader.headerOk())
    {
        controlPointsPtr_.set
        (
            new vectorIOField
            (
                IOobject
                (
                    "controlPoints",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            )
        );
    }
    else
    {
        controlPointsPtr_.set
        (
            new vectorIOField
            (
                IOobject
                (
                    "controlPoints",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                aMesh().areaCentres().internalField()
            )
        );

        initializeControlPointsPosition();
    }
}

void Foam::interfaceTrackingFvMesh::makeTotalDisplacement()
{
    if (debug)
    {
        Info<< "surfaceTracking::makeTotalDisplacement() : "
            << "making zero total points displacement"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!totalDisplacementPtr_.empty())
    {
        FatalErrorIn("surfaceTracking::makeTotalDisplacement()")
            << "total points displacement already exists"
            << abort(FatalError);
    }

    totalDisplacementPtr_.set
    (
        new vectorIOField
        (
            IOobject
            (
                "totalDisplacement",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            vectorField
            (
                mesh().boundaryMesh()[surfacePatchID()].nPoints(),
                vector::zero
            )
        )
    );
}

void Foam::interfaceTrackingFvMesh::makeMotionPointsMask()
{
    DebugInFunction
        << "making motion points mask" << nl;

    if (!motionPointsMaskPtr_.empty())
    {
        FatalErrorInFunction
            << "motion points mask already exists"
            << abort(FatalError);
    }

    motionPointsMaskPtr_.set
    (
        new scalarField
        (
            mesh().boundaryMesh()[surfacePatchID()].nPoints(),
            1
        )
    );

    // Mark free surface boundary points
    // that belong to processor patches
    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == processorFaPatch::typeName
        )
        {
            const labelList& patchPoints =
                aMesh().boundary()[patchI].pointLabels();

            forAll(patchPoints, pointI)
            {
                motionPointsMask()[patchPoints[pointI]] = -1;
            }
        }
    }

    // Mark fixed free surface boundary points
    for (const word& patchName : fixedSurfacePatches_)
    {
        const label fixedPatchID = aMesh().boundary().findPatchID(patchName);

        if (fixedPatchID == -1)
        {
            FatalErrorInFunction
                << "Wrong faPatch name in the fixedSurfacePatches list"
                << " defined in the dynamicMeshDict dictionary"
                << abort(FatalError);
        }

        const labelList& patchPoints =
            aMesh().boundary()[fixedPatchID].pointLabels();

        forAll(patchPoints, pointI)
        {
            motionPointsMask()[patchPoints[pointI]] = 0;
        }
    }
}


void Foam::interfaceTrackingFvMesh::makeDirections()
{
    DebugInFunction
        << "make displacement directions for points and control points" << nl;

    if (!pointsDisplacementDirPtr_.empty() || !facesDisplacementDirPtr_.empty())
    {
        FatalErrorInFunction
            << "points, control points displacement directions already exist"
            << abort(FatalError);
    }

    pointsDisplacementDirPtr_.set
    (
        new vectorField
        (
            mesh().boundaryMesh()[surfacePatchID()].nPoints(),
            pTraits<vector>::zero
        )
    );

    facesDisplacementDirPtr_.set
    (
        new vectorField
        (
            mesh().boundaryMesh()[surfacePatchID()].size(),
            pTraits<vector>::zero
        )
    );

    if (!normalMotionDir())
    {
        if (mag(motionDir_) < SMALL)
        {
            FatalErrorInFunction
                << "Zero motion direction"
                << abort(FatalError);
        }

        facesDisplacementDir() = motionDir_;
        pointsDisplacementDir() = motionDir_;
    }

    updateDisplacementDirections();
}


void Foam::interfaceTrackingFvMesh::updateDisplacementDirections()
{
    if(normalMotionDir())
    {
        // Update point displacement correction
        pointsDisplacementDir() = aMesh().pointAreaNormals();

        // Correct point displacement direction
        // at the "centerline" symmetryPlane which represents the axis
        // of an axisymmetric case
        forAll(aMesh().boundary(), patchI)
        {
            if(aMesh().boundary()[patchI].type() == wedgeFaPatch::typeName)
            {
                const wedgeFaPatch& wedgePatch =
                    refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

                vector axis = wedgePatch.axis();

                label centerLinePatchID =
                    aMesh().boundary().findPatchID("centerline");

                if(centerLinePatchID != -1)
                {
                    const labelList& pointLabels =
                        aMesh().boundary()[centerLinePatchID].pointLabels();

                    forAll(pointLabels, pointI)
                    {
                        vector dir =
                            pointsDisplacementDir()[pointLabels[pointI]];

                        dir = (dir&axis)*axis;
                        dir /= mag(dir);

                        pointsDisplacementDir()[pointLabels[pointI]] = dir;
                    }
                }
                else
                {
                    Info << "Warning: centerline polyPatch does not exist. "
                        << "Surface points displacement directions "
                        << "will not be corrected at the axis (centerline)"
                        << endl;
                }

//                break;
            }
        }

        label nWedgePatches = 0;
        vector wedgeDirVec = vector::zero;
        forAll(mesh().boundaryMesh(), patchI)
        {
            if (isA<wedgePolyPatch>(mesh().boundaryMesh()[patchI]))
            {
                const wedgePolyPatch& wpp = refCast<const wedgePolyPatch>
                (
                    mesh().boundaryMesh()[patchI]
                );

                nWedgePatches++;
                wedgeDirVec += cmptMag(wpp.centreNormal());
            }
        }

        reduce(nWedgePatches, maxOp<label>());

        if (nWedgePatches)
        {
            reduce(wedgeDirVec, sumOp<vector>());

            wedgeDirVec /= mag(wedgeDirVec);

            vectorField dir = pointsDisplacementDir();

            for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
            {
                if (wedgeDirVec[cmpt] > 1e-6)
                {
                    dir.replace(cmpt, 0.0);
                }
            }

            dir /= mag(dir);

            pointsDisplacementDir() = dir;
        }

        // project the point displacement dirs to the surface on which the A phase
        // is moving if this surface is a wall patch
//        cLine().projectPointDisplDirToWall();

        // Update face displacement direction
        facesDisplacementDir() =
            aMesh().faceAreaNormals().internalField();

        // Correction of control points postion
        const vectorField& Cf = aMesh().areaCentres().internalField();

        controlPoints() =
            facesDisplacementDir()
          * (facesDisplacementDir()&(controlPoints() - Cf))
          + Cf;
    }
}


void Foam::interfaceTrackingFvMesh::initializeControlPointsPosition()
{
    scalarField deltaH = scalarField(aMesh().nFaces(), 0.0);

    pointField displacement = pointDisplacement(deltaH);

    const faceList& faces = aMesh().faces();
    const pointField& points = aMesh().points();

    pointField newPoints = points + displacement;

    scalarField sweptVol(faces.size(), 0.0);

    forAll(faces, faceI)
    {
        sweptVol[faceI] = -faces[faceI].sweptVol(points, newPoints);
    }

    vectorField faceArea(faces.size(), vector::zero);

    forAll (faceArea, faceI)
    {
        faceArea[faceI] = faces[faceI].normal(newPoints);
    }

    forAll(deltaH, faceI)
    {
        deltaH[faceI] = sweptVol[faceI]/
            (faceArea[faceI] & facesDisplacementDir()[faceI]);
    }

    forAll(fixedSurfacePatches_, patchI)
    {
        label fixedPatchID = 
            aMesh().boundary().findPatchID
            (
                fixedSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("interfaceTrackingFvMesh::initializeControlPointsPosition()")
                << "Wrong faPatch name in the fixedSurfacePatches list"
                    << " defined in the dynamicMeshDict dictionary"
                    << abort(FatalError);
        }

        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        forAll(eFaces, edgeI)
        {
            deltaH[eFaces[edgeI]] *= 2.0;
        }
    }

    displacement = pointDisplacement(deltaH);
}


void Foam::interfaceTrackingFvMesh::correctPointDisplacement
(
    const scalarField& sweptVolCorr,
    vectorField& displacement
)
{
    const labelListList& pFaces =
        aMesh().patch().pointFaces();

    const faceList& faces =
        aMesh().patch().localFaces();

    const pointField& points =
        aMesh().patch().localPoints();

    forAll(fixedSurfacePatches_, patchI)
    {
        label fixedPatchID =
            aMesh().boundary().findPatchID
            (
                fixedSurfacePatches_[patchI]
            );

        const labelList& pLabels =
            aMesh().boundary()[fixedPatchID].pointLabels();

        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        labelHashSet pointSet;

        forAll(eFaces, edgeI)
        {
            label curFace = eFaces[edgeI];

            const labelList& curPoints = faces[curFace];
            
            forAll(curPoints, pointI)
            {
                label curPoint = curPoints[pointI];
                label index = findIndex(pLabels, curPoint);

                if (index == -1)
                {
                    if (!pointSet.found(curPoint))
                    {
                        pointSet.insert(curPoint);
                    }
                }
            }
        }

        labelList corrPoints = pointSet.toc();

        labelListList corrPointFaces(corrPoints.size());

        forAll(corrPoints, pointI)
        {
            label curPoint = corrPoints[pointI];

            labelHashSet faceSet;
            
            forAll(pFaces[curPoint], faceI)
            {
                label curFace = pFaces[curPoint][faceI];

                label index = findIndex(eFaces, curFace);

                if (index != -1)
                {
                    if (!faceSet.found(curFace))
                    {
                        faceSet.insert(curFace);
                    }
                }
            }

            corrPointFaces[pointI] = faceSet.toc();
        }

        forAll(corrPoints, pointI)
        {
            label curPoint = corrPoints[pointI];

            scalar curDisp = 0;

            const labelList& curPointFaces = corrPointFaces[pointI];

            forAll(curPointFaces, faceI)
            {
                const face& curFace = faces[curPointFaces[faceI]];

                label ptInFace = curFace.which(curPoint);
                label next = curFace.nextLabel(ptInFace);
                label prev = curFace.prevLabel(ptInFace);
                
                vector a = points[next] - points[curPoint];
                vector b = points[prev] - points[curPoint];
                const vector& c = pointsDisplacementDir()[curPoint];

                curDisp += 2*sweptVolCorr[curPointFaces[faceI]]/((a^b)&c);
            }

            curDisp /= curPointFaces.size();

            displacement[curPoint] = 
                curDisp*pointsDisplacementDir()[curPoint];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceTrackingFvMesh::interfaceTrackingFvMesh
(
    const IOobject& io
)
:
//    dynamicMotionSolverFvMesh(io),
    topoChangerFvMesh(io),
    motionPtr_(motionSolver::New(*this)),
    aMeshPtr_(new faMesh(*this)),
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
    movingSurfacePatches_(),
    surfacePatchID_ //(-1),
    (
        this->boundaryMesh().findPatchID("freeSurface")
    ),
    fixedSurfacePatches_(),
    nonReflectingFreeSurfacePatches_(),
    pointNormalsCorrectionPatches_(),
    normalMotionDir_
    (
        motionDict_.lookup("normalMotionDir")
    ),
    motionDir_(pTraits<vector>::zero),
    smoothing_(false),
    rigidFreeSurface_
    (
        motionDict_.lookupOrDefault<Switch>("rigidFreeSurface", false)
    ),
    timeIndex_(-1),
    sweptVolCorrOld_(this->boundaryMesh()[surfacePatchID_].size(), 0.0),
    resetFluxFrequency_(100),
    controlPointsPtr_(),
    totalDisplacementPtr_(),
    motionPointsMaskPtr_(),
    pointsDisplacementDirPtr_(),
    facesDisplacementDirPtr_()
//    contactAnglePtr_(nullptr)
{
    initializeData();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceTrackingFvMesh::~interfaceTrackingFvMesh()
{}


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


Foam::vectorField& Foam::interfaceTrackingFvMesh::controlPoints()
{
    if (controlPointsPtr_.empty())
    {
        makeControlPoints();
    }

    return controlPointsPtr_();
}

Foam::vectorField& Foam::interfaceTrackingFvMesh::totalDisplacement()
{
    if (totalDisplacementPtr_.empty())
    {
        makeTotalDisplacement();
    }

    return totalDisplacementPtr_();
}

Foam::scalarField& Foam::interfaceTrackingFvMesh::motionPointsMask()
{
    if (motionPointsMaskPtr_.empty())
    {
        makeMotionPointsMask();
    }

    return motionPointsMaskPtr_();
}


Foam::vectorField& Foam::interfaceTrackingFvMesh::pointsDisplacementDir()
{
    if (pointsDisplacementDirPtr_.empty())
    {
        makeDirections();
    }

    return pointsDisplacementDirPtr_();
}


Foam::vectorField& Foam::interfaceTrackingFvMesh::facesDisplacementDir()
{
    if (facesDisplacementDirPtr_.empty())
    {
        makeDirections();
    }

    return facesDisplacementDirPtr_();
}


bool Foam::interfaceTrackingFvMesh::update()
{
    if (timeIndex_ != mesh().time().timeIndex())
    {
        if (smoothing_ && !rigidFreeSurface_)
        {
            clearControlPoints();
        }

        updateDisplacementDirections();

        const scalarField& K = aMesh().faceCurvatures().internalField();

        Info<< "Free surface curvature: min = " << gMin(K)
            << ", max = " << gMax(K) << ", average = " << gAverage(K) << nl;

        timeIndex_ = mesh().time().timeIndex();
    }

    if (!rigidFreeSurface_)
    {
        // This is currently relative flux
        scalarField sweptVolCorr =
            phi().boundaryField()[surfacePatchID()];

//        if (mesh().moving())
//        {
//            sweptVolCorr -=
//                fvc::meshPhi(U())().boundaryField()[surfacePatchID()];
//        }

        pointField newMeshPoints = mesh().allPoints();

        Info<< "Moving surface continuity error : sum local = "
            << gSum(mag(sweptVolCorr)) << ", global = " << gSum(sweptVolCorr)
            << endl;

        word ddtScheme
        (
            mesh().schemesDict().ddtScheme
            (
                "ddt(" + U().name() + ')'
            )
        );

        if
        (
            ddtScheme
         == fv::CrankNicolsonDdtScheme<vector>::typeName
        )
        {
            sweptVolCorr *= (1.0/2.0)*mesh().time().deltaT().value();
        }
        else if
        (
            ddtScheme
         == fv::EulerDdtScheme<vector>::typeName
        )
        {
            sweptVolCorr *= mesh().time().deltaT().value();
        }
        else if
        (
            (ddtScheme == "bdf2")
         ||
            (ddtScheme == fv::backwardDdtScheme<vector>::typeName)
        )
        {
            if (mesh().time().timeIndex() == 1)
            {
                sweptVolCorr *= mesh().time().deltaT().value();
            }
            else
            {
                sweptVolCorr *= (2.0/3.0)*mesh().time().deltaT().value();

//                sweptVolCorr += (1.0/3.0)*sweptVolCorrOld_;
            }
        }
        else
        {
            FatalErrorInFunction
                << "Unsupported temporal differencing scheme : "
                << ddtScheme << nl
                << abort(FatalError);
        }

        const scalarField& Sf = aMesh().S();
        const vectorField& Nf = aMesh().faceAreaNormals().internalField();

        scalarField deltaHf
        (
            sweptVolCorr/(Sf*(Nf & facesDisplacementDir()))
        );

        forAll(fixedSurfacePatches_, patchI)
        {
            label fixedPatchID = 
                aMesh().boundary().findPatchID
                (
                    fixedSurfacePatches_[patchI]
                );

            if(fixedPatchID == -1)
            {
                FatalErrorIn("interfaceTrackingFvMesh::update(...)")
                    << "Wrong faPatch name in the fixedSurfacePatches list"
                        << " defined in the surfaceProperties dictionary"
                        << abort(FatalError);
            }
            
            const labelList& eFaces =
                aMesh().boundary()[fixedPatchID].edgeFaces();

            forAll(eFaces, edgeI)
            {
                deltaHf[eFaces[edgeI]] *= 2.0;
            }
        }

        controlPoints() += facesDisplacementDir()*deltaHf;
//        updateDisplacementDirections();

        pointField displacement(pointDisplacement(deltaHf)); //line 291 in freeSurface

        //- correct point displacement for fixed surface patches
        correctPointDisplacement(sweptVolCorr, displacement);

        // Move only free surface points
        const labelList& meshPoints =
            mesh().boundaryMesh()[surfacePatchID()].meshPoints();

        forAll (displacement, pointI)
        {
            newMeshPoints[meshPoints[pointI]] += displacement[pointI];
        }

        // Update total displacement field
//        totalDisplacement() += displacement;

        if (timeIndex_ < mesh().time().timeIndex())
        {
            totalDisplacement() = displacement;

            timeIndex_ = mesh().time().timeIndex();
        }
        else
        {
            totalDisplacement() += displacement;
        }

        // mesh motion
//        fvMotionSolver& mSolver =
//            dynamic_cast<fvMotionSolver&>
//            (
//                motionPtr_()
//            );

//        pointVectorField& motionU = mSolver.pointMotionU();

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
                motionU.boundaryField()[surfacePatchID()]
            );

//        motionUPatch ==
//            totalDisplacement()/mesh().time().deltaT().value();

        motionUPatch ==
            displacement/mesh().time().deltaT().value();


//        twoDPointCorrector twoDPointCorr(mesh());
//        twoDPointCorr.correctPoints(newMeshPoints);

        sweptVolCorrOld_ = sweptVolCorr;

//        fvMesh::movePoints(newMeshPoints);
//        setOldPoints(newMeshPoints);
//        movePoints(newMeshPoints);

//        aMesh().movePoints();

        motionPtr_->solve();
//-        mSolver.solve();

//        dynamicMotionSolverFvMesh::update();

        fvMesh::movePoints(motionPtr_->curPoints());
//-        fvMesh::movePoints(mSolver.curPoints());

//        fvMesh::movePoints(mSolver.newPoints());
    }
    else
    {
        notImplemented
        (
            "interfaceTrackingFvMesh::update()\n"
            "Rigid free-surface model\n"
            "not implemented"
        );
    }

    // Enforce interface motion based on actual flux (without artefacts),
    // by reset of residual fluxes across free surface
    if
    (
        ((timeIndex_-1) % resetFluxFrequency_) == 0
     && timeIndex_ > 1
    )
    {
        //- Non-const access to flux on patch
        fvsPatchField<scalar>& surfacePhiField = 
            const_cast<fvsPatchField<scalar>& >
            (
                mesh().lookupObject<surfaceScalarField>("phi")
                .boundaryField()[surfacePatchID()]
            );

        surfacePhiField *= 0;
    }

    // dynamicMotionSolverFvMesh::update();

    return true;
}


void Foam::interfaceTrackingFvMesh::writeVTK() const
{
    // Write patch and points into VTK
    OFstream mps(mesh().time().timePath()/"surface.vtk");

    const vectorField& points = aMesh().patch().points();
    const IndirectList<face>& faces = aMesh().patch();

    mps << "# vtk DataFile Version 2.0" << nl
        << mesh().time().timePath()/"surface.vtk" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << points.size() << " float" << nl;

    // Write points
    List<float> mlpBuffer(3*points.size());

    label counter = 0;
    forAll(points, i)
    {
        mlpBuffer[counter++] = float(points[i].x());
        mlpBuffer[counter++] = float(points[i].y());
        mlpBuffer[counter++] = float(points[i].z());
    }

    forAll(mlpBuffer, i)
    {
        mps << mlpBuffer[i] << ' ';

        if (i > 0 && (i % 10) == 0)
        {
            mps << nl;
        }
    }

    // Write faces
    label nFaceVerts = 0;

    forAll(faces, faceI)
    {
        nFaceVerts += faces[faceI].size() + 1;
    }
    labelList mlfBuffer(nFaceVerts);

    counter = 0;
    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        mlfBuffer[counter++] = f.size();

        forAll(f, fpI)
        {
            mlfBuffer[counter++] = f[fpI];
        }
    }
    mps << nl;

    mps << "POLYGONS " << faces.size() << ' ' << nFaceVerts << endl;

    forAll(mlfBuffer, i)
    {
        mps << mlfBuffer[i] << ' ';

        if (i > 0 && (i % 10) == 0)
        {
            mps << nl;
        }
    }
    mps << nl;

     aMesh().patch().writeVTK
     (
         mesh().time().timePath()/"surface",
         aMesh().patch(),
         aMesh().patch().points()
     );
}


void Foam::interfaceTrackingFvMesh::writeVTKControlPoints()
{
    // Write control points into VTK
    fileName name(mesh().time().timePath()/"surfaceControlPoints.vtk");
    OFstream mps(name);

    Info<< "Writing surface control point to " << name << endl;

    mps << "# vtk DataFile Version 2.0" << nl
        << name << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << controlPoints().size() << " float" << nl;

    forAll(controlPoints(), pointI)
    {
        mps << controlPoints()[pointI].x() << ' '
            << controlPoints()[pointI].y() << ' '
            << controlPoints()[pointI].z() << nl;
    }

    // Write vertices
    mps << "VERTICES " << controlPoints().size() << ' '
        << controlPoints().size()*2 << nl;

    forAll(controlPoints(), pointI)
    {
        mps << 1 << ' ' << pointI << nl;
    }
}


// ************************************************************************* //
