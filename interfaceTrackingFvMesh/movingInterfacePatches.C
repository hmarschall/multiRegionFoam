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

\*---------------------------------------------------------------------------*/

#include "movingInterfacePatches.H"
#include "polyTopoChanger.H"
#include "foamTime.H"
#include "addToRunTimeSelectionTable.H"

#include "polyPatchID.H"
#include "wedgeFaPatch.H"
#include "wedgeFaPatchFields.H"
#include "slipFaPatchFields.H"
#include "fixedValueFaPatchFields.H"
#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"
#include "wallFvPatch.H"
#include "fvcMeshPhi.H"

#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"

#include "coordinateSystem.H"
#include "scalarMatrices.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(movingInterfacePatches, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::movingInterfacePatches::initializeData()
{
    if (motionDict_.found("neighbourMesh"))
    {
        isInterface_ = true;
    }

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
        const word surfacePatchName(name_);

        polyPatchID patch(surfacePatchName, mesh_.boundaryMesh());

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


void Foam::movingInterfacePatches::makeControlPoints()
{
    if (debug)
    {
        Info<< "movingInterfacePatches::makeControlPoints() : "
            << "making control points"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!controlPointsPtr_.empty())
    {
        FatalErrorIn("movingInterfacePatches::makeInterpolators()")
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

void Foam::movingInterfacePatches::makeTotalDisplacement()
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
                mesh().boundaryMesh()[patchID()].nPoints(),
                vector::zero
            )
        )
    );
}

void Foam::movingInterfacePatches::makeMotionPointsMask()
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
            mesh().boundaryMesh()[patchID()].nPoints(),
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


void Foam::movingInterfacePatches::makeDirections()
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
            mesh().boundaryMesh()[patchID()].nPoints(),
            pTraits<vector>::zero
        )
    );

    facesDisplacementDirPtr_.set
    (
        new vectorField
        (
            mesh().boundaryMesh()[patchID()].size(),
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


void Foam::movingInterfacePatches::updateDisplacementDirections()
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

        // TODO: adapt to new framework
        // project the point displacement dirs to the surface on which 
        // the A phase is moving if this surface is a wall patch
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


void Foam::movingInterfacePatches::initializeControlPointsPosition()
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
            FatalErrorIn("movingInterfacePatches::initializeControlPointsPosition()")
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


void Foam::movingInterfacePatches::correctPointDisplacement
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


void Foam::movingInterfacePatches::makeGlobalPatches() const
{
    if (globalPatchPtr_.valid() || globalNbrPatchPtr_.valid())
    {
        FatalErrorIn(type() + "::makeGlobalPatches() const")
            << "Pointer already set!" << abort(FatalError);
    }

    Info<< "Creating global patches : "
        << patch().name() << "and "
        << nbrPatch().name()
        << endl;

    globalPatchPtr_.set(new globalPolyPatch(patch().name(), mesh()));
    globalNbrPatchPtr_.set(new globalPolyPatch(nbrPatch().name(), nbrMesh()));
}

void Foam::movingInterfacePatches::clearGlobalPatches() const
{
    globalPatchPtr_.clear();
    globalNbrPatchPtr_.clear();
}


void Foam::movingInterfacePatches::makeInterfaceToInterface() const
{
    if (interfaceToInterfacePtr_.valid())
    {
        FatalErrorIn
        (
            "void Foam::movingInterfacePatches::"
            "makeInterfaceToInterface() const"
        )   << "Mapping object already set!" << abort(FatalError);
    }

    // Lookup the type 
    const word type = motionDict_.lookupOrDefault<word>
    (
        "interfaceTransferMethod", "GGI"
    );

    interfaceToInterfacePtr_ =
    (
        interfaceToInterfaceMapping::New
        (
            type,
            motionDict_.subDict(type + "Coeffs"),
            mesh().boundaryMesh()[patchID()],
            nbrMesh().boundaryMesh()[nbrPatch().index()],
            globalPatch(),
            globalNbrPatch()
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingInterfacePatches::movingInterfacePatches
(
    const word& name,
    const dynamicFvMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    aMeshPtr_(new faMesh(mesh)),
    motionDict_(dict),
    surfacePatchID_
    (
        mesh_.boundaryMesh().findPatchID(name_)
    ),
    surfacePatch_(mesh_.boundary()[surfacePatchID_]),
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
    isInterface_(false),
    motionTrigger_(true),
    interfaceDeformationLimit_
    (
        motionDict_.lookupOrDefault<dimensionedScalar>
        (
            "interfaceDeformationLimit",
            dimensionedScalar("interfaceDeformationLimit", dimLength, 0.0)
        )
    ),
    timeIndex_(-1),
    sweptVolCorrOld_(mesh_.boundaryMesh()[surfacePatchID_].size(), 0.0),
    resetFluxFrequency_(100),
    neighbourRegionName_
    (
        motionDict_
        .lookupOrDefault<word>("neighbourMesh", word::null)
    ),
    neighbourPatchName_
    (
        motionDict_
        .lookupOrDefault<word>("neighbourPatch", word::null)
    ),
    controlPointsPtr_(),
    totalDisplacementPtr_(),
    motionPointsMaskPtr_(),
    pointsDisplacementDirPtr_(),
    facesDisplacementDirPtr_(),
//    contactAnglePtr_(nullptr)
    globalPatchPtr_(),
    globalNbrPatchPtr_(),
    interfaceToInterfacePtr_(),
    interpolatorUpdateFrequency_
    (
        motionDict_
        .lookupOrDefault<int>("interpolatorUpdateFrequency", 0)
    )
{
    initializeData();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::movingInterfacePatches::~movingInterfacePatches()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingInterfacePatches::updateTopology()
{
    controlPointsPtr_.clear();
    motionPointsMaskPtr_.clear();
    pointsDisplacementDirPtr_.clear();
    facesDisplacementDirPtr_.clear();
    totalDisplacementPtr_.clear();

    if (isInterface_)
    {
        globalPatchPtr_().updateMesh();
        globalNbrPatchPtr_().updateMesh();
    }
}

const Foam::dynamicFvMesh& Foam::movingInterfacePatches::nbrMesh() const
{
    return
    (
        mesh_.objectRegistry::parent()
        .lookupObject<dynamicFvMesh>(neighbourRegionName_)
    );
}

const Foam::fvPatch& Foam::movingInterfacePatches::nbrPatch() const
{
    label neighbourPatchID = 
        nbrMesh().boundaryMesh().findPatchID(neighbourPatchName_);

    return
    (
        nbrMesh().boundary()[neighbourPatchID]
    );
}

Foam::vectorField& Foam::movingInterfacePatches::controlPoints()
{
    if (controlPointsPtr_.empty())
    {
        makeControlPoints();
    }

    return controlPointsPtr_();
}

Foam::vectorField& Foam::movingInterfacePatches::totalDisplacement()
{
    if (totalDisplacementPtr_.empty())
    {
        makeTotalDisplacement();
    }

    return totalDisplacementPtr_();
}

Foam::scalarField& Foam::movingInterfacePatches::motionPointsMask()
{
    if (motionPointsMaskPtr_.empty())
    {
        makeMotionPointsMask();
    }

    return motionPointsMaskPtr_();
}


Foam::vectorField& Foam::movingInterfacePatches::pointsDisplacementDir()
{
    if (pointsDisplacementDirPtr_.empty())
    {
        makeDirections();
    }

    return pointsDisplacementDirPtr_();
}


Foam::vectorField& Foam::movingInterfacePatches::facesDisplacementDir()
{
    if (facesDisplacementDirPtr_.empty())
    {
        makeDirections();
    }

    return facesDisplacementDirPtr_();
}


Foam::tmp<vectorField> 
Foam::movingInterfacePatches::surfacePointDisplacement()
{
    if (timeIndex_ == -1) //initial only / do not move into constructor
    {
        if (motionDict_.found("neighbourMesh"))
        {
            // Create global patches synced in parallel runs such that
            // all faces are present on all processors
            makeGlobalPatches();

            // Force creation of interface-to-interface object 
            // as they may need to read fields on restart
            interfaceToInterface();
        }
    }

    const pointField& points = aMesh().patch().localPoints();

    tmp<vectorField> tdisplacement
    (
        new vectorField
        (
            points.size(),
            vector::zero
        )
    );

    vectorField& displacement = tdisplacement();

    if (timeIndex_ != mesh().time().timeIndex())
    {
        if (smoothing_ && !rigidFreeSurface_)
        {
            clearControlPoints();
        }

        updateDisplacementDirections();

        if (isInterface_)
        {
            updateInterpolatorAndGlobalPatches();
        }

        const scalarField& K = aMesh().faceCurvatures().internalField();

        Info<< "Free surface curvature: min = " << gMin(K)
            << ", max = " << gMax(K) << ", average = " << gAverage(K) << nl;

        timeIndex_ = mesh().time().timeIndex();
    }

    if (!rigidFreeSurface_)
    {
        const surfaceScalarField& phi =
            mesh().lookupObject<surfaceScalarField>("phi");

        // This is currently relative flux
        scalarField sweptVolCorr =
            phi.boundaryField()[patchID()];

//        if (mesh().moving())
//        {
//            sweptVolCorr -=
//                fvc::meshPhi(U())().boundaryField()[patchID()];
//        }

        pointField newMeshPoints = mesh().allPoints();

        Info<< "Moving surface continuity error : sum local = "
            << gSum(mag(sweptVolCorr)) << ", global = " << gSum(sweptVolCorr)
            << endl;

        const volVectorField& U = 
            mesh().lookupObject<volVectorField>("U");

        word ddtScheme
        (
            mesh().schemesDict().ddtScheme
            (
                "ddt(" + U.name() + ')'
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

                //sweptVolCorr += (1.0/3.0)*sweptVolCorrOld_;
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

        displacement = pointDisplacement(deltaHf); //line 291 in freeSurface

        //- correct point displacement for fixed surface patches
        correctPointDisplacement(sweptVolCorr, displacement);

        sweptVolCorrOld_ = sweptVolCorr;

        // TODO: can newMeshPoints and totalDisplacement be removed?
        // Move only free surface points
        const labelList& meshPoints =
            mesh().boundaryMesh()[patchID()].meshPoints();

        forAll (displacement, pointI)
        {
            newMeshPoints[meshPoints[pointI]] += displacement[pointI];
        }

        // Update total displacement field
        //totalDisplacement() += displacement;

        if (timeIndex_ < mesh().time().timeIndex())
        {
            totalDisplacement() = displacement;

            timeIndex_ = mesh().time().timeIndex();
        }
        else
        {
            totalDisplacement() += displacement;
        }

        // Move whole mesh only if interface deformation limit is exceeded
        // Note: This check is supposed to be done for every interface in 
        // the list, thus triggering mesh motion automatically if the 
        // individual limit is exceeded by ONE of the moving interfaces.
        scalar minCellThickness =
            2*gMin(1.0/mesh().boundary()[patchID()].deltaCoeffs());

        scalar maxInterfaceDeformation =
            gMax(mag(totalDisplacement()))/minCellThickness;

        Info << "Maximal relative interface deformation: "
            << maxInterfaceDeformation << endl;

        if 
        (
            !(maxInterfaceDeformation > interfaceDeformationLimit_.value())
        )
        {
            motionTrigger_ = false;

            return tmp<vectorField>
            (
                new vectorField
                (
                    points.size(),
                    vector::zero
                )
            );
        }
    }
    else
    {
        notImplemented
        (
            "movingInterfacePatches::surfacePointDisplacement()\n"
            "Rigid free-surface model\n"
            "not implemented"
        );
    }


    // Enforce interface motion based on actual flux (without artefacts),
    // by reset of residual fluxes across free surface
//    if
//    (
//        ((timeIndex_-1) % resetFluxFrequency_) == 0
//     && timeIndex_ > 1
//    )
//    {
//        //- Non-const access to flux on patch
//        fvsPatchField<scalar>& surfacePhiField = 
//            const_cast<fvsPatchField<scalar>& >
//            (
//                mesh().lookupObject<surfaceScalarField>("phi")
//                .boundaryField()[patchID()]
//            );

//        surfacePhiField *= 0;
//    }

    return tdisplacement;
}

Foam::tmp<vectorField> 
Foam::movingInterfacePatches::shadowPointDisplacement
(
    const vectorField& displacement
)
{
    tmp<vectorField> shadowDisplacement
    (
        new vectorField(displacement.size())
    );

    // Interpolate displacement to global patch point data
    vectorField globalDisplacement =
        globalPatch().patchPointToGlobal(displacement);

    // Create global displacement field of same size
    vectorField globalShadowDisplacement(globalDisplacement.size());

    // Interpolate globally patch-to-patch
    interfaceToInterface().transferPointsZoneToZone
    (
        globalPatch().globalPatch(),        // from zone
        globalNbrPatch().globalPatch(),     // to zone
        globalDisplacement,                 // from field
        globalShadowDisplacement            // to field
    );

    // Move global patch points due to mesh motion
    globalPatchPtr_().movePoints(globalDisplacement);
    globalNbrPatchPtr_().movePoints(globalShadowDisplacement);

    // Filter global patch point data to patch
    shadowDisplacement() = globalNbrPatch().globalPointToPatch
    (
        globalShadowDisplacement
    );

    // Return
    return shadowDisplacement;
}

void movingInterfacePatches::correctPointNormals()
{
    // Correct normals for fixed patches points
    Info << "Correct point normals" << endl;

    vectorField& N =
        const_cast<vectorField&>
        (
            aMesh().pointAreaNormals()
        );

    const labelListList& pFaces =
        aMesh().patch().pointFaces();

    const labelListList& fFaces =
        aMesh().patch().faceFaces();

    const faceList& faces = 
        aMesh().patch().localFaces();

    const pointField& points = 
        aMesh().patch().localPoints();

    // Wedge points
//    forAll(aMesh().boundary(), patchI)
//    {
//        if (aMesh().boundary()[patchI].type() == wedgeFaPatch::typeName)
//        {
//            const wedgeFaPatch& wedgePatch =
//                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

//            const labelList& patchPoints = wedgePatch.pointLabels();

//            forAll(patchPoints, pointI)
//            {
//                label curPoint = patchPoints[pointI];

//                labelHashSet faceSet;
//                forAll(pFaces[curPoint], faceI)
//                {
//                    faceSet.insert(pFaces[curPoint][faceI]);
//                }
//                labelList curFaces = faceSet.toc();
//                
//                labelHashSet pointSet;

//                pointSet.insert(curPoint);
//                for(label i=0; i<curFaces.size(); i++)
//                {
//                    const labelList& facePoints = faces[curFaces[i]];
//                    for(label j=0; j<facePoints.size(); j++)
//                    {
//                        if(!pointSet.found(facePoints[j]))
//                        {
//                            pointSet.insert(facePoints[j]);
//                        }
//                    }
//                }
//                pointSet.erase(curPoint);
//                labelList curPoints = pointSet.toc();

//                
//                labelHashSet addPointsSet;
//                forAll(curPoints, pointI)
//                {
//                    label index = 
//                        findIndex(patchPoints, curPoints[pointI]);

//                    if (index != -1)
//                    {
//                        addPointsSet.insert(curPoints[pointI]);
//                    }  
//                }
//                addPointsSet.insert(curPoint);
//                labelList curAddPoints = addPointsSet.toc();


//                if (curPoints.size() + curAddPoints.size() >= 5)
//                {
//                    vectorField allPoints
//                    (
//                        curPoints.size()+curAddPoints.size()
//                    );
//                    scalarField W(curPoints.size()+curAddPoints.size(), 1.0);
//                    label I = -1;
//                    for(label i=0; i<curPoints.size(); i++)
//                    {
//                        I++;
//                        allPoints[I] = points[curPoints[i]];
//                        W[I] = 1.0/magSqr(allPoints[I] - points[curPoint]);
//                    }
//                    for(label i=0; i<curAddPoints.size(); i++)
//                    {
//                        I++;
//                        allPoints[I] = 
//                            transform
//                            (
//                                wedgePatch.faceT(),
//                                points[curAddPoints[i]]
//                            );
//                        W[I] = 1.0/magSqr(allPoints[I] - points[curPoint]);
//                    }

//                    // Transforme points
//                    vector origin = points[curPoint];
//                    vector axis = N[curPoint]/mag(N[curPoint]);
//                    vector dir = (allPoints[0] - points[curPoint]);
//                    dir -= axis*(axis&dir);
//                    dir /= mag(dir);
//                    coordinateSystem cs("cs", origin, axis, dir);
//                    
//                    forAll(allPoints, pI)
//                    {
//                        allPoints[pI] = cs.localPosition(allPoints[pI]);
//                    }
//                    
//                    scalarRectangularMatrix M
//                    (
//                        allPoints.size(),
//                        5,
//                        0.0
//                    );

//                    for(label i = 0; i < allPoints.size(); i++)
//                    {
//                        M[i][0] = sqr(allPoints[i].x());
//                        M[i][1] = sqr(allPoints[i].y());
//                        M[i][2] = allPoints[i].x()*allPoints[i].y();
//                        M[i][3] = allPoints[i].x();
//                        M[i][4] = allPoints[i].y();
//                    }
//                    
//                    scalarSquareMatrix MtM(5, 0.0);

//                    for (label i = 0; i < MtM.n(); i++)
//                    {
//                        for (label j = 0; j < MtM.m(); j++)
//                        {
//                            for (label k = 0; k < M.n(); k++)
//                            {
//                                MtM[i][j] += M[k][i]*M[k][j]*W[k];
//                            }
//                        }
//                    }
//                    
//                    scalarField MtR(5, 0);

//                    for (label i=0; i<MtR.size(); i++)
//                    {
//                        for (label j=0; j<M.n(); j++)
//                        {
//                            MtR[i] += M[j][i]*allPoints[j].z()*W[j];
//                        }
//                    }
//            
//                    scalarSquareMatrix::LUsolve(MtM, MtR);

//                    vector curNormal = vector(MtR[3], MtR[4], -1);
//                    
//                    curNormal = cs.globalVector(curNormal);
//                    
//                    curNormal *= sign(curNormal&N[curPoint]);
//                    
//                    N[curPoint] = curNormal;
//                }
//            }
//        }
//    }
    
    // Fixed boundary points

    forAll(fixedSurfacePatches_, patchI)
    {
        label fixedPatchID = 
            aMesh().boundary().findPatchID
            (
                fixedSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("movingInterfacePatches::correctPointNormals()")
                << "Wrong faPatch name in the fixedSurfacePatches list"
                    << " defined in the surfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& pLabels = 
            aMesh().boundary()[fixedPatchID].pointLabels();

        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();
 
        forAll(pLabels, pointI)
        {
            label curPoint = pLabels[pointI];

            labelHashSet faceSet;
            forAll(pFaces[curPoint], faceI)
            {
                faceSet.insert(pFaces[curPoint][faceI]);
            }

            labelList curFaces = faceSet.toc();

            forAll(curFaces, faceI)
            {
                const labelList& curFaceFaces =
                    fFaces[curFaces[faceI]];
                
                forAll(curFaceFaces, fI)
                {
                    label curFaceFace = curFaceFaces[fI];
                        
                    label index = findIndex(eFaces, curFaceFace);

                    if( (index==-1) && !faceSet.found(curFaceFace) )
                    {
                        faceSet.insert(curFaceFace);
                    }
                }
            }
            curFaces = faceSet.toc();

            labelHashSet pointSet;

            pointSet.insert(curPoint);
            for(label i=0; i<curFaces.size(); i++)
            {
                const labelList& fPoints = faces[curFaces[i]];
                for(label j=0; j<fPoints.size(); j++)
                {
                    if(!pointSet.found(fPoints[j]))
                    {
                        pointSet.insert(fPoints[j]);
                    }
                }
            }

            pointSet.erase(curPoint);

            labelList curPoints = pointSet.toc();

            // LS quadric fit
            vectorField allPoints(curPoints.size());
            scalarField W(curPoints.size(), 1.0);
            for(label i=0; i<curPoints.size(); i++)
            {
                allPoints[i] = points[curPoints[i]];
                W[i] = 1.0/magSqr(allPoints[i] - points[curPoint]);
            }

            // Transforme points
            vector origin = points[curPoint];
            vector axis = N[curPoint]/mag(N[curPoint]);
            vector dir = (allPoints[0] - points[curPoint]);
            dir -= axis*(axis&dir);
            dir /= mag(dir);
            coordinateSystem cs("cs", origin, axis, dir);

            forAll(allPoints, pI)
            {
                allPoints[pI] = cs.localPosition(allPoints[pI]);
            }

            scalarRectangularMatrix M
            (
                allPoints.size(),
                5,
                0.0
            );

            for(label i = 0; i < allPoints.size(); i++)
            {
                M[i][0] = sqr(allPoints[i].x());
                M[i][1] = sqr(allPoints[i].y());
                M[i][2] = allPoints[i].x()*allPoints[i].y();
                M[i][3] = allPoints[i].x();
                M[i][4] = allPoints[i].y();
            }

            scalarSquareMatrix MtM(5, 0.0);

            for (label i = 0; i < MtM.n(); i++)
            {
                for (label j = 0; j < MtM.m(); j++)
                {
                    for (label k = 0; k < M.n(); k++)
                    {
                        MtM[i][j] += M[k][i]*M[k][j]*W[k];
                    }
                }
            }

            scalarField MtR(5, 0);

            for (label i=0; i<MtR.size(); i++)
            {
                for (label j=0; j<M.n(); j++)
                {
                    MtR[i] += M[j][i]*allPoints[j].z()*W[j];
                }
            }

            scalarSquareMatrix::LUsolve(MtM, MtR);

            vector curNormal = vector(MtR[3], MtR[4], -1);

            curNormal = cs.globalVector(curNormal);

            curNormal *= sign(curNormal&N[curPoint]);

            N[curPoint] = curNormal;
        }
    }

    // Correcte wedge points
//    forAll (aMesh().boundary(), patchI)
//    {
//        if (aMesh().boundary()[patchI].type() == wedgeFaPatch::typeName)
//        {
//            const wedgeFaPatch& wedgePatch =
//                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

//            const labelList& patchPoints = wedgePatch.pointLabels();

//            vector n =
//                transform
//                (
//                    wedgePatch.edgeT(),
//                    wedgePatch.centreNormal()
//                );

//            n /= mag(n);

//            forAll (patchPoints, pointI)
//            {
//                N[patchPoints[pointI]]
//                    -= n*(n&N[patchPoints[pointI]]);
//            }
//        }
//    }


    // Boundary points correction
    forAll (aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().correctPatchPointNormals(patchI) 
        && !aMesh().boundary()[patchI].coupled()
        )
        {
            if (aMesh().boundary()[patchI].ngbPolyPatchIndex() == -1)
            {
                FatalErrorIn
                    (
                        "void correctPointNormals::correctPointNormals()"
                    )   << "Neighbour polyPatch index is not defined "
                        << "for faPatch " << aMesh().boundary()[patchI].name()
                        << abort(FatalError);
            }

            labelList patchPoints = aMesh().boundary()[patchI].pointLabels();

            vectorField n = 
                aMesh().boundary()[patchI].ngbPolyPatchPointNormals();

            forAll (patchPoints, pointI)
            {
                N[patchPoints[pointI]]
                    -= n[pointI]*(n[pointI]&N[patchPoints[pointI]]);
            }
        }
    }


    N /= mag(N);
}

const Foam::interfaceToInterfaceMapping&
Foam::movingInterfacePatches::interfaceToInterface() const
{
    if (interfaceToInterfacePtr_.empty())
    {
        makeInterfaceToInterface();
    }

    return interfaceToInterfacePtr_();
}


const Foam::globalPolyPatch& Foam::movingInterfacePatches::globalPatch() const
{
    if (globalPatchPtr_.empty())
    {
        FatalErrorIn(type() + "::makeGlobalPatch() const")
            << "makeGlobalPatches(...) must be called before globalPatchA "
            << "can be called!" << abort(FatalError);
    }

    return globalPatchPtr_();
}

const Foam::globalPolyPatch& Foam::movingInterfacePatches::globalNbrPatch() const
{
    if (globalNbrPatchPtr_.empty())
    {
        FatalErrorIn(type() + "::makeGlobalPatch() const")
            << "makeGlobalPatches(...) must be called before globalPatchA "
            << "can be called!" << abort(FatalError);
    }

    return globalNbrPatchPtr_();
}


void Foam::movingInterfacePatches::updateInterpolatorAndGlobalPatches()
{
    Info << "Updating interpolator and global patches" << endl;

    if (interfaceToInterfacePtr_.empty())
    {
        interfaceToInterface();
    }
    else if (interpolatorUpdateFrequency_ != 0)
    {
        if
        (
            (
                (mesh().time().timeIndex() - 1) % interpolatorUpdateFrequency_
            ) == 0
        )
        {
            // Re-create interpolators
            interfaceToInterface();

            // Enforce topology update
            updateTopology();
        }
    }
}


void Foam::movingInterfacePatches::writeVTK() const
{
    // Write patch and points into VTK
    OFstream mps(mesh().time().timePath()/name_+".vtk");

    const vectorField& points = aMesh().patch().points();
    const IndirectList<face>& faces = aMesh().patch();

    mps << "# vtk DataFile Version 2.0" << nl
        << mesh().time().timePath()/name_+".vtk" << nl
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


void Foam::movingInterfacePatches::writeVTKControlPoints()
{
    // Write control points into VTK
    fileName name(mesh().time().timePath()/name_+"ControlPoints.vtk");
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
