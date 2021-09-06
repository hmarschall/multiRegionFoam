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
        dynamicFvMesh,
        interfaceTrackingFvMesh,
        IOobject
    );
//    addToRunTimeSelectionTable
//    (
//        dynamicFvMesh,
//        interfaceTrackingFvMesh,
//        doInit
//    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//void Foam::interfaceTrackingFvMesh::initializeData()
//{
//    // Set free surface patch index
//    {
//        const word fsPatchName(motion().get<word>("fsPatchName"));

//        polyPatchID patch(fsPatchName, this->boundaryMesh());

//        if (!patch.active())
//        {
//            FatalErrorInFunction
//                << "Patch name " << fsPatchName << " not found."
//                << abort(FatalError);
//        }

//        surfacePatchID_ = patch.index();
//    }

//    // Set point normal correction for finite area mesh
//    {
//        boolList& correction = aMesh().correctPatchPointNormals();

//        for (const word& patchName : pointNormalsCorrectionPatches_)
//        {
//            label patchID = aMesh().boundary().findPatchID(patchName);

//            if (patchID == -1)
//            {
//                FatalErrorInFunction
//                    << "Patch name '" << patchName
//                    << "' for point normals correction does not exist"
//                    << abort(FatalError);
//            }

//            correction[patchID] = true;
//        }
//    }

//    // Read motion direction
//    if (!normalMotionDir_)
//    {
//        motionDir_ = normalised(motion().get<vector>("motionDir"));
//    }

//    // Check if contact angle is defined
//    makeContactAngle();

//    motion().readIfPresent
//    (
//        "nonReflectingFreeSurfacePatches",
//        nonReflectingFreeSurfacePatches_
//    );
//}


void Foam::interfaceTrackingFvMesh::makeSurfaceNetPhi() const
{
    DebugInFunction
        << "making surface net flux" << nl;

    if (surfaceNetPhiPtr_)
    {
        FatalErrorInFunction
            << "surface net flux already exists"
            << abort(FatalError);
    }

    surfaceNetPhiPtr_ = new areaScalarField
    (
        IOobject
        (
            "surfaceNetPhi",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh(),
        dimensionedScalar("surfaceNetPhi0", dimVelocity*dimArea, pTraits<scalar>::zero)
    );
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
    if (controlPointsPtr_)
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
        controlPointsPtr_ =
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
            );
    }
    else
    {
        controlPointsPtr_ =
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
            );

        initializeControlPointsPosition();
    }
}


void Foam::interfaceTrackingFvMesh::makeMotionPointsMask()
{
    DebugInFunction
        << "making motion points mask" << nl;

    if (motionPointsMaskPtr_)
    {
        FatalErrorInFunction
            << "motion points mask already exists"
            << abort(FatalError);
    }

    motionPointsMaskPtr_ = new scalarField//labelList
    (
        mesh().boundaryMesh()[surfacePatchID()].nPoints(),
        1
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
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
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

    if (pointsDisplacementDirPtr_ || facesDisplacementDirPtr_)
    {
        FatalErrorInFunction
            << "points, control points displacement directions already exist"
            << abort(FatalError);
    }

    pointsDisplacementDirPtr_ =
        new vectorField
        (
            mesh().boundaryMesh()[surfacePatchID()].nPoints(),
            pTraits<vector>::zero
        );

    facesDisplacementDirPtr_ =
        new vectorField
        (
            mesh().boundaryMesh()[surfacePatchID()].size(),
            pTraits<vector>::zero
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

//    forAll(fixedSurfacePatches_, patchI)
//    {
//        label fixedPatchID = 
//            aMesh().boundary().findPatchID
//            (
//                fixedSurfacePatches_[patchI]
//            );

//        if(fixedPatchID == -1)
//        {
//            FatalErrorIn("surfaceTracking::surfaceTracking(...)")
//                << "Wrong faPatch name in the fixedSurfacePatches list"
//                    << " defined in the surfaceProperties dictionary"
//                    << abort(FatalError);
//        }

//        const labelList& eFaces =
//            aMesh().boundary()[fixedPatchID].edgeFaces();

//        forAll(eFaces, edgeI)
//        {
//            deltaH[eFaces[edgeI]] *= 2.0;
//        }
//    }

    displacement = pointDisplacement(deltaH);
}


//Foam::scalar Foam::interfaceTrackingFvMesh::maxCourantNumber()
//{
//    scalar CoNum = 0;

//    if (pureFreeSurface())
//    {
//        const scalarField& dE = aMesh().lPN();

//        CoNum = gMax
//        (
//            mesh().time().deltaT().value()/
//            sqrt
//            (
//                Foam::pow(dE, 3.0)/2.0/M_PI/(sigma().value() + SMALL)
//            )
//        );
//    }
//    else
//    {
//        scalarField sigmaE
//        (
//            linearEdgeInterpolate(surfaceTension())().internalField().field()
//          + SMALL
//        );

//        const scalarField& dE = aMesh().lPN();

//        CoNum = gMax
//        (
//            mesh().time().deltaT().value()/
//            sqrt
//            (
//                Foam::pow(dE, 3.0)/2.0/M_PI/sigmaE
//            )
//        );
//    }

//    return CoNum;
//}


void Foam::interfaceTrackingFvMesh::updateProperties()
{
    const singlePhaseTransportModel& properties =
        mesh().lookupObject<singlePhaseTransportModel>
        (
            "transportProperties"
        );

    rho_ = dimensionedScalar(properties.lookup("rho"));

    // volScalarField mu_ = rho_*properties.nu();

    sigma0_ = dimensionedScalar(properties.lookup("sigma"));
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
    dynamicMotionSolverFvMesh(io),
    aMeshPtr_(new faMesh(*this)),
    surfacePatchID_(-1),
    fixedSurfacePatches_(),
    nonReflectingFreeSurfacePatches_(),
    pointNormalsCorrectionPatches_(),
    normalMotionDir_(false),
    motionDir_(pTraits<vector>::zero),
    smoothing_(false),
    pureFreeSurface_(true),
    rigidFreeSurface_(false),
    correctContactLineNormals_(false),
    sigma0_("zero", dimForce/dimLength/dimDensity, pTraits<scalar>::zero),
    rho_("one", dimDensity, 1.0),
    timeIndex_(-1),
    controlPointsPtr_(nullptr),
    motionPointsMaskPtr_(nullptr),
    pointsDisplacementDirPtr_(nullptr),
    facesDisplacementDirPtr_(nullptr),
    surfaceNetPhiPtr_(nullptr)
//    surfaceTensionPtr_(nullptr),
//    contactAnglePtr_(nullptr)
{
//    if (doInit)
//    {
//        init(false);    // do not initialise lower levels
//    }
}

/*
Foam::interfaceTrackingFvMesh::interfaceTrackingFvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    labelList&& allOwner,
    labelList&& allNeighbour,
    const bool syncPar
)
:
    dynamicMotionSolverFvMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(allOwner),
        std::move(allNeighbour),
        syncPar
    ),
    aMeshPtr_(new faMesh(*this)),
    surfacePatchID_(-1),
    fixedSurfacePatches_(),
    nonReflectingFreeSurfacePatches_(),
    pointNormalsCorrectionPatches_(),
    normalMotionDir_(false),
    motionDir_(pTraits<vector>::zero),
    smoothing_(false),
    pureFreeSurface_(true),
    sigma0_("zero", dimForce/dimLength/dimDensity, pTraits<scalar>::zero),
    rho_("one", dimDensity, 1.0),
    timeIndex_(-1),
    controlPointsPtr_(nullptr),
    motionPointsMaskPtr_(nullptr),
    pointsDisplacementDirPtr_(nullptr),
    facesDisplacementDirPtr_(nullptr),
    surfaceNetPhiPtr_(nullptr),
    surfaceTensionPtr_(nullptr),
    contactAnglePtr_(nullptr)
{}
*/

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceTrackingFvMesh::~interfaceTrackingFvMesh()
{
    deleteDemandDrivenData(controlPointsPtr_);
    deleteDemandDrivenData(motionPointsMaskPtr_);
    deleteDemandDrivenData(pointsDisplacementDirPtr_);
    deleteDemandDrivenData(facesDisplacementDirPtr_);
    deleteDemandDrivenData(surfaceNetPhiPtr_);
//    deleteDemandDrivenData(surfaceTensionPtr_);
//    deleteDemandDrivenData(contactAnglePtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//bool Foam::interfaceTrackingFvMesh::init(const bool doInit)
//{
//    if (doInit)
//    {
//        dynamicMotionSolverFvMesh::init(doInit);
//    }

//    aMeshPtr_.reset(new faMesh(*this));

//    // Set motion-based data
//    fixedSurfacePatches_ =
//        motion().get<wordList>("fixedFreeSurfacePatches");

//    pointNormalsCorrectionPatches_ =
//        motion().get<wordList>("pointNormalsCorrectionPatches");

//    normalMotionDir_ = motion().get<bool>("normalMotionDir");
//    smoothing_ = motion().getOrDefault("smoothing", false);
//    pureFreeSurface_ = motion().getOrDefault("pureFreeSurface", true);

//    initializeData();

//    return true;
//}


Foam::areaScalarField& Foam::interfaceTrackingFvMesh::surfaceNetPhi()
{
    if (!surfaceNetPhiPtr_)
    {
        makeSurfaceNetPhi();
    }

    return *surfaceNetPhiPtr_;
}


const Foam::areaScalarField& Foam::interfaceTrackingFvMesh::surfaceNetPhi() const
{
    if (!surfaceNetPhiPtr_)
    {
        makeSurfaceNetPhi();
    }

    return *surfaceNetPhiPtr_;
}


const Foam::volVectorField& Foam::interfaceTrackingFvMesh::U() const
{
    return mesh().objectRegistry::lookupObject<volVectorField>("U");
}


//const Foam::volScalarField& Foam::interfaceTrackingFvMesh::p() const
//{
//    return *this.objectRegistry::lookupObject<volScalarField>("p");
//}


const Foam::surfaceScalarField& Foam::interfaceTrackingFvMesh::phi() const
{
    return mesh().objectRegistry::lookupObject<surfaceScalarField>("phi");
}


Foam::vectorField& Foam::interfaceTrackingFvMesh::controlPoints()
{
    if (!controlPointsPtr_)
    {
        makeControlPoints();
    }

    return *controlPointsPtr_;
}


//Foam::labelList& Foam::interfaceTrackingFvMesh::motionPointsMask()
//{
//    if (!motionPointsMaskPtr_)
//    {
//        makeMotionPointsMask();
//    }

//    return *motionPointsMaskPtr_;
//}


Foam::vectorField& Foam::interfaceTrackingFvMesh::pointsDisplacementDir()
{
    if (!pointsDisplacementDirPtr_)
    {
        makeDirections();
    }

    return *pointsDisplacementDirPtr_;
}


Foam::vectorField& Foam::interfaceTrackingFvMesh::facesDisplacementDir()
{
    if (!facesDisplacementDirPtr_)
    {
        makeDirections();
    }

    return *facesDisplacementDirPtr_;
}


//Foam::areaScalarField&
//Foam::interfaceTrackingFvMesh::surfaceTension()
//{
//    if (!surfaceTensionPtr_)
//    {
//        makeSurfaceTension();
//    }

//    return *surfaceTensionPtr_;
//}


//const Foam::areaScalarField&
//Foam::interfaceTrackingFvMesh::surfaceTension() const
//{
//    if (!surfaceTensionPtr_)
//    {
//        makeSurfaceTension();
//    }

//    return *surfaceTensionPtr_;
//}


bool Foam::interfaceTrackingFvMesh::update()
{
    if (timeIndex_ != mesh().time().timeIndex())
    {
        if (smoothing_ && !rigidFreeSurface_)
        {
            clearControlPoints();
        }

        updateDisplacementDirections();

        updateProperties();

        Info<< "Maximal capillary Courant number: "
            << maxCourantNumber() << endl;

        const scalarField& K = aMesh().faceCurvatures().internalField();

        Info<< "Free surface curvature: min = " << gMin(K)
            << ", max = " << gMax(K) << ", average = " << gAverage(K) << nl;

        timeIndex_ = mesh().time().timeIndex();
    }

//    if (!rigidFreeSurface_)
//    {
        pointField newMeshPoints = mesh().allPoints();

        // This is currently relative flux
        scalarField sweptVolCorr =
            phi().boundaryField()[surfacePatchID()];

        // Info<< "Free surface flux: sum local = "
        //     << gSum(mag(sweptVolCorr))
        //     << ", global = " << gSum(sweptVolCorr) << endl;

        // if (mesh().moving())
        // {
        //     sweptVolCorr -=
        //         fvc::meshPhi(U())().boundaryField()[fsPatchIndex()];
        // }

        Info<< "Free surface continuity error : sum local = "
            << gSum(mag(sweptVolCorr)) << ", global = " << gSum(sweptVolCorr)
            << endl;

        // For postprocessing
        surfaceNetPhi().internalField() = sweptVolCorr;

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
            ddtScheme
         == fv::backwardDdtScheme<vector>::typeName
        )
        {
            if (mesh().time().timeIndex() == 1)
            {
                sweptVolCorr *= mesh().time().deltaT().value();
            }
            else
            {
                sweptVolCorr *= (2.0/3.0)*mesh().time().deltaT().value();
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
                    << "Wrong faPatch name in the fixedFreeSurfacePatches list"
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

        pointField displacement(pointDisplacement(deltaHf));
        correctPointDisplacement(sweptVolCorr, displacement);

        //- HM
        // Move only free surface points
        const labelList& meshPoints =
            mesh().boundaryMesh()[surfacePatchID()].meshPoints();

        forAll (displacement, pointI)
        {
            newMeshPoints[meshPoints[pointI]] += displacement[pointI];
        }

        // TODO: add totalDisplacementPtr_ data
        // OR can this be replace by displacement? (see line 178 in freeSurface.C)
        // totalDisplacement() += displacement;

        twoDPointCorrector twoDPointCorr(mesh());

        twoDPointCorr.correctPoints(newMeshPoints);

        mesh().movePoints(newMeshPoints);


        // mesh motion
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

        motionUPatch ==
            displacement/mesh().time().deltaT().value();

//        velocityMotionSolver& vMotion =
//            refCast<velocityMotionSolver>
//            (
//                const_cast<motionSolver&>(motion())
//            );

//        pointVectorField& pointMotionU = vMotion.pointMotionU();
//        pointMotionU.primitiveFieldRef() = pTraits<vector>::zero;

//        fixedValuePointPatchVectorField& fsPatchPointMeshU =
//            refCast<fixedValuePointPatchVectorField>
//            (
//                const_cast<pointPatchVectorField&>
//                (
//                    pointMotionU.boundaryField()[fsPatchIndex()]
//                )
//            );

//        fsPatchPointMeshU ==
//            displacement/mesh().time().deltaT().value();

        dynamicMotionSolverFvMesh::update();

//    }
//    else
//    {
//        vectorField displacement
//        (
//            mesh().boundaryMesh()[fsPatchIndex()].nPoints(),
//            pTraits<vector>::zero
//        );

//        velocityMotionSolver& vMotion =
//            refCast<velocityMotionSolver>
//            (
//                const_cast<motionSolver&>(motion())
//            );

//        pointVectorField& pointMotionU = vMotion.pointMotionU();
//        pointMotionU.primitiveFieldRef() = pTraits<vector>::zero;

//        fixedValuePointPatchVectorField& fsPatchPointMeshU =
//            refCast<fixedValuePointPatchVectorField>
//            (
//                const_cast<pointPatchVectorField&>
//                (
//                    pointMotionU.boundaryField()[fsPatchIndex()]
//                )
//            );

//        fsPatchPointMeshU ==
//            displacement/mesh().time().deltaT().value();

//        dynamicMotionSolverFvMesh::update();
//    }

    return true;
}


//void Foam::interfaceTrackingFvMesh::writeVTK() const
//{
//    // Write patch and points into VTK
//    OFstream mps(mesh().time().timePath()/"freeSurface.vtk");

//    const vectorField& points = aMesh().patch().points();
//    const IndirectList<face>& faces = aMesh().patch();

//    mps << "# vtk DataFile Version 2.0" << nl
//        << mesh().time().timePath()/"freeSurface.vtk" << nl
//        << "ASCII" << nl
//        << "DATASET POLYDATA" << nl
//        << "POINTS " << points.size() << " float" << nl;

//    // Write points
//    List<float> mlpBuffer(3*points.size());

//    label counter = 0;
//    forAll(points, i)
//    {
//        mlpBuffer[counter++] = float(points[i].x());
//        mlpBuffer[counter++] = float(points[i].y());
//        mlpBuffer[counter++] = float(points[i].z());
//    }

//    forAll(mlpBuffer, i)
//    {
//        mps << mlpBuffer[i] << ' ';

//        if (i > 0 && (i % 10) == 0)
//        {
//            mps << nl;
//        }
//    }

//    // Write faces
//    label nFaceVerts = 0;

//    forAll(faces, faceI)
//    {
//        nFaceVerts += faces[faceI].size() + 1;
//    }
//    labelList mlfBuffer(nFaceVerts);

//    counter = 0;
//    forAll(faces, faceI)
//    {
//        const face& f = faces[faceI];

//        mlfBuffer[counter++] = f.size();

//        forAll(f, fpI)
//        {
//            mlfBuffer[counter++] = f[fpI];
//        }
//    }
//    mps << nl;

//    mps << "POLYGONS " << faces.size() << ' ' << nFaceVerts << endl;

//    forAll(mlfBuffer, i)
//    {
//        mps << mlfBuffer[i] << ' ';

//        if (i > 0 && (i % 10) == 0)
//        {
//            mps << nl;
//        }
//    }
//    mps << nl;

//    // aMesh().patch().writeVTK
//    // (
//    //     mesh().time().timePath()/"freeSurface",
//    //     aMesh().patch(),
//    //     aMesh().patch().points()
//    // );
//}


//void Foam::interfaceTrackingFvMesh::writeVTKControlPoints()
//{
//    // Write control points into VTK
//    fileName name(mesh().time().timePath()/"freeSurfaceControlPoints.vtk");
//    OFstream mps(name);

//    Info<< "Writing free surface control point to " << name << endl;

//    mps << "# vtk DataFile Version 2.0" << nl
//        << name << nl
//        << "ASCII" << nl
//        << "DATASET POLYDATA" << nl
//        << "POINTS " << controlPoints().size() << " float" << nl;

//    forAll(controlPoints(), pointI)
//    {
//        mps << controlPoints()[pointI].x() << ' '
//            << controlPoints()[pointI].y() << ' '
//            << controlPoints()[pointI].z() << nl;
//    }

//    // Write vertices
//    mps << "VERTICES " << controlPoints().size() << ' '
//        << controlPoints().size()*2 << nl;

//    forAll(controlPoints(), pointI)
//    {
//        mps << 1 << ' ' << pointI << nl;
//    }
//}


// ************************************************************************* //
