/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Zeljko Tukovic, FSB Zagreb.
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

#include "movingInterfacePatches.H"
#include "primitivePatchInterpolation.H"
#include "emptyFaPatch.H"
#include "wedgeFaPatch.H"
#include "wallFvPatch.H"
#include "PstreamCombineReduceOps.H"
#include "coordinateSystem.H"
#include "unitConversion.H"
#include "scalarMatrices.H"
#include "tensor2D.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//Foam::tmp<Foam::vectorField>
//Foam::movingInterfacePatches::pointDisplacement(const scalarField& deltaH)
//{
//    const pointField& points = aMesh().patch().localPoints();
//    const labelListList& pointFaces = aMesh().patch().pointFaces();

//    const labelList faceCells = 
//        mesh().boundary()[patchID()].patch().faceCells();

//    controlPoints() += facesDisplacementDir()*deltaH;

//    // Correct for curvature at axis
////    if (Pstream::master())
//    {
//        label patchID = aMesh().boundary().findPatchID("centerline");

//        if (patchID != -1)
//        {
//            const labelList& eFaces =
//                aMesh().boundary()[patchID].edgeFaces();

//            const labelListList& fFaces = aMesh().patch().faceFaces();
//            const vectorField& fCentres =
//                aMesh().areaCentres().internalField();

//            forAll(eFaces, edgeI)
//            {
//                const label& curFace = eFaces[edgeI];
//                const labelList& curFaceFaces = fFaces[curFace];

//                scalar H = 0.0;
//                label counter = 0;

//                forAll(curFaceFaces, faceI)
//                {
//                    label index = findIndex(eFaces, curFaceFaces[faceI]);

//                    if (index == -1)
//                    {
//                        H +=
//                            facesDisplacementDir()[curFaceFaces[faceI]]
//                          & (
//                                controlPoints()[curFaceFaces[faceI]]
//                              - fCentres[curFaceFaces[faceI]]
//                            );

//                        counter++;
//                    }
//                }

//                H /= counter;

//                controlPoints()[curFace] =
//                    fCentres[curFace]
//                  + facesDisplacementDir()[curFace]*H;
//            }
//        }
//    }

//    // Correct controPoints next to fixed patches
//    if (Pstream::master())
//    {
//        forAll(fixedSurfacePatches_, patchI)
//        {
//            label fixedPatchID =
//                aMesh().boundary().findPatchID
//                (
//                    fixedSurfacePatches_[patchI]
//                );

//            if(fixedPatchID == -1)
//            {
//                FatalErrorIn("freeSurface::freeSurface(...)")
//                    << "Wrong faPatch name in the fixedSurfacePatches list"
//                        << " defined in the freeSurfaceProperties dictionary"
//                        << abort(FatalError);
//            }

//            const labelList& eFaces =
//                aMesh().boundary()[fixedPatchID].edgeFaces();

//            const labelListList& fFaces = aMesh().patch().faceFaces();
//            const vectorField& fCentres =
//                aMesh().areaCentres().internalField();
//            
//            forAll(eFaces, edgeI)
//            {
//                const label& curFace = eFaces[edgeI];
//                const labelList& curFaceFaces = fFaces[curFace];

//                scalar H = 0.0;
//                label counter = 0;

//                forAll(curFaceFaces, faceI)
//                {
//                    label index = findIndex(eFaces, curFaceFaces[faceI]);
//                    
//                    if (index == -1)
//                    {
//                        H +=
//                            facesDisplacementDir()[curFaceFaces[faceI]]
//                          & (
//                                controlPoints()[curFaceFaces[faceI]]
//                              - fCentres[curFaceFaces[faceI]]
//                            );

//                        counter++;
//                    }
//                }
//                
//                H /= counter;

//                controlPoints()[curFace] =
//                    fCentres[curFace]
//                  + facesDisplacementDir()[curFace]*H;
//            }
//        }
//    }

//    // Calculate displacement of internal points
//    tmp<vectorField> tdisplacement
//    (
//        new vectorField
//        (
//            points.size(),
//            vector::zero
//        )
//    );

//    vectorField& displacement = tdisplacement();

//    forAll (pointFaces, pointI)
//    {
//        scalar weightsSum = 0.0;
//        const labelList& curPointFaces = pointFaces[pointI];

//        forAll (curPointFaces, faceI)
//        {
//            label curFace = curPointFaces[faceI];

//            scalar weight = 1.0/mag
//            (
//                points[pointI]
//              - controlPoints()[curFace]
//            );

//            displacement[pointI] += weight*controlPoints()[curFace];

//            weightsSum += weight;
//        }

//        displacement[pointI] /= weightsSum;

//        displacement[pointI] -= points[pointI];
//    }

//    displacement = motionPointsMask()*
//        (pointsDisplacementDir()&displacement)*
//        pointsDisplacementDir();

//    // Calculate displacement of axis point
//    forAll (aMesh().boundary(), patchI)
//    {
//        if
//        (
//            aMesh().boundary()[patchI].type()
//         == wedgeFaPatch::typeName
//        )
//        {
//            const wedgeFaPatch& wedgePatch =
//                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

////            if(wedgePatch.axisPoint() > -1)
//            forAll(wedgePatch.axisPoints(), apI)
//            {
////                label axisPoint = wedgePatch.axisPoint();
//                label axisPoint = wedgePatch.axisPoints()[apI];

//                displacement[axisPoint] =
//                    pointsDisplacementDir()[axisPoint]
//                   *(
//                        pointsDisplacementDir()[axisPoint]
//                       &(
//                            controlPoints()[pointFaces[axisPoint][0]]
//                          - points[axisPoint]
//                        )
//                    );
//            }
//        }
//    }

//    return tdisplacement;
//}


//tmp<vectorField> movingInterfacePatches::pointDisplacement(const scalarField& deltaH) 
//{
//    const pointField& points = aMesh().patch().localPoints();
//    const labelListList& pointFaces = aMesh().patch().pointFaces();

//    controlPoints() += facesDisplacementDir()*deltaH;

//    tmp<vectorField> tdisplacement
//    (
//        new vectorField
//        (
//            points.size(),
//            vector::zero
//        )
//    );

//    vectorField& displacement = tdisplacement();


//    // Calculate displacement of internal points

//    const edgeList& edges = aMesh().patch().edges();

//    labelList internalPoints = aMesh().internalPoints();

//    forAll (internalPoints, pointI)
//    {
//        label curPoint = internalPoints[pointI];

//        const labelList& curPointFaces = pointFaces[curPoint];

//        scalarField w(curPointFaces.size(), 0.0);

//        forAll (curPointFaces, faceI)
//        {
//            label curFace = curPointFaces[faceI];

//            scalar magDistance =
//            (
//                 stabilise
//                 (
//                     mag
//                     (
//                        controlPoints()[curFace]
//                      - points[curPoint]
//                     ), VSMALL
//                 )
//            );

//            w[faceI] = 1.0/magDistance;
//        }

//        w /= sum(w);

//        vector Q = vector::zero;

//        forAll (curPointFaces, faceI)
//        {
//            label curFace = curPointFaces[faceI];

//            Q += w[faceI]*controlPoints()[curFace];
//        }

//        displacement[curPoint] =
//        (
//           -1.0*
//           (
//               pointsDisplacementDir()[curPoint] & (points[curPoint] - Q)
//              *pointsDisplacementDir()[curPoint]
//           )
//        );
//    }

//    // Calculate displacement of points
//    // which belonge to empty and wedge patches
//    labelList boundaryPoints = aMesh().boundaryPoints();

//    const labelListList& edgeFaces = aMesh().patch().edgeFaces();
//    const labelListList& pointEdges = aMesh().patch().pointEdges();

//    const vectorField& pointNormals = aMesh().pointAreaNormals();

//    forAll (boundaryPoints, pointI)
//    {
//        label curPoint = boundaryPoints[pointI];

//        if (motionPointsMask()[curPoint])
//        {
//            // Calculating mirror points
//            const labelList& curPointEdges = pointEdges[curPoint];

//            vectorField mirrorPoints(2, vector::zero);

//            label counter = -1;

//            forAll (curPointEdges, edgeI)
//            {
//                label curEdge = curPointEdges[edgeI];

//                if(edgeFaces[curEdge].size() == 1)
//                {
//                    label patchID = -1;

//                    forAll(aMesh().boundary(), patchI)
//                    {
//                        // Note: aMesh().boundary()[patchI].size()==0 for
//                        // empty faPatch
//                        forAll(aMesh().boundary()[patchI], eI)
//                        {
//                            if (aMesh().boundary()[patchI][eI] == curEdge)
//                            {
//                                patchID = patchI;
//                                break;
//                            }
//                        }
//                    }

//                    if
//                    (
//                        patchID > -1
//                     && (
//                            aMesh().boundary()[patchID].type()
//                         == wedgeFaPatch::typeName
//                        )
//                    )
//                    {
//                        const wedgeFaPatch& wedgePatch =
//                            refCast<const wedgeFaPatch>
//                            (
//                                aMesh().boundary()[patchID]
//                            );

//                        mirrorPoints[++counter] =
//                            transform
//                            (
//                                wedgePatch.faceT(),
//                                controlPoints()[edgeFaces[curEdge][0]]
//                            );
//                    }
//                    else
//                    {
//                        vector nE =
//                            pointNormals[edges[curEdge].start()]
//                          + pointNormals[edges[curEdge].end()];

//                        nE /= mag(nE);

//                        vector eP =
//                            controlPoints()[edgeFaces[curEdge][0]]
//                          - edges[curEdge].centre(points);

//                        mirrorPoints[++counter] =
//                            edges[curEdge].centre(points)
//                          - ((I - 2.0*nE*nE)&eP);
//                    }
//                }
//            }


//            // Calculating LS plane interpolation
//            const labelList& curPointFaces = pointFaces[curPoint];

//            symmTensor M = symmTensor::zero;

//            vector S = vector::zero;

//            scalarField w(curPointFaces.size() + 2, 0.0);

//            forAll (curPointFaces, faceI)
//            {
//                label curFace = curPointFaces[faceI];

//                w[faceI] = 1.0/mag
//                (
//                    controlPoints()[curFace]
//                  - points[curPoint]
//                );
//            }

//            forAll (mirrorPoints, pI)
//            {
//                w[curPointFaces.size() + pI] = 1.0/mag
//                (
//                    mirrorPoints[pI]
//                  - points[curPoint]
//                );
//            }

//            w /= sum(w);


//            forAll (curPointFaces, faceI)
//            {
//                label curFace = curPointFaces[faceI];

//                M = M + sqr(w[faceI])*sqr(controlPoints()[curFace]);

//                S += sqr(w[faceI])*controlPoints()[curFace];
//            }


//            forAll (mirrorPoints, pI)
//            {
//                M = M + sqr(w[curPointFaces.size()+pI])*sqr(mirrorPoints[pI]);

//                S += sqr(w[curPointFaces.size()+pI])*mirrorPoints[pI];
//            }

//            vector N = inv(M)&S;

//            N /= mag(N);

//            scalar p = (S&N)/sum(sqr(w));

//            displacement[curPoint] =
//                pointsDisplacementDir()[curPoint]*
//                (p - (points[curPoint]&N))/
//                (pointsDisplacementDir()[curPoint]&N);
//        }
//    }

//    forAll(aMesh().boundary(), patchI)
//    {
//        bool fixedPatch = false;

//        forAll(fixedSurfacePatches_, fpI)
//        {
//            label fixedPatchID = aMesh().boundary().findPatchID
//            (
//                fixedSurfacePatches_[fpI]
//            );

//            if(fixedPatchID == -1)
//            {
//                FatalErrorIn("surfaceTracking::pointDisplacemen(...)")
//                    << "Wrong faPatch name in the fixedSurfacePatches list"
//                        << " defined in the surfaceProperties dictionary"
//                        << abort(FatalError);
//            }

//            if (fixedPatchID == patchI)
//            {
//                fixedPatch = true;
//            }
//        }

//        if
//        (
//            (
//                aMesh().boundary()[patchI].type()
//             != emptyFaPatch::typeName
//            )
//         && (
//                aMesh().boundary()[patchI].type()
//             != wedgeFaPatch::typeName
//            )
//         && (
//                aMesh().boundary()[patchI].type()
//             != processorFaPatch::typeName
//            )
//         && !fixedPatch
//        )
//        {
//            labelList patchPoints =
//                aMesh().boundary()[patchI].pointLabels();

//            labelListList patchPointEdges =
//                aMesh().boundary()[patchI].pointEdges();

//            unallocLabelList patchEdgeFaces =
//                aMesh().boundary()[patchI].edgeFaces();

//            forAll(patchPoints, pointI)
//            {
//                forAll(patchPointEdges[pointI], edgeI)
//                {
//                    label curEdge = patchPointEdges[pointI][edgeI];

//                    displacement[patchPoints[pointI]] +=
//                        pointsDisplacementDir()[patchPoints[pointI]]
//                       *(
//                            pointsDisplacementDir()[patchPoints[pointI]]
//                           &(
//                                controlPoints()[patchEdgeFaces[curEdge]]
//                              - points[patchPoints[pointI]]
//                            )
//                        );
//                }

//                displacement[patchPoints[pointI]] /=
//                    patchPointEdges[pointI].size();
//            }
//        }
//    }


//    // Calculate displacement of axis point
//    forAll (aMesh().boundary(), patchI)
//    {
//        if
//        (
//            aMesh().boundary()[patchI].type()
//         == wedgeFaPatch::typeName
//        )
//        {
//            const wedgeFaPatch& wedgePatch =
//                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

//            if(wedgePatch.axisPoint() > -1)
//            {
//                label axisPoint = wedgePatch.axisPoint();

//                displacement[axisPoint] =
//                    pointsDisplacementDir()[axisPoint]
//                   *(
//                        pointsDisplacementDir()[axisPoint]
//                       &(
//                            controlPoints()[pointFaces[axisPoint][0]]
//                          - points[axisPoint]
//                        )
//                    );
//            }
//        }
//    }


//    // Calculate displacement of processor patch points
//    forAll (aMesh().boundary(), patchI)
//    {
//        if
//        (
//            aMesh().boundary()[patchI].type()
//         == processorFaPatch::typeName
//        )
//        {
//            const processorFaPatch& procPatch =
//                refCast<const processorFaPatch>(aMesh().boundary()[patchI]);

//            const labelList& patchPointLabels =
//                procPatch.pointLabels();

//            Field<symmTensor> ownM(patchPointLabels.size(), symmTensor::zero);
//            vectorField ownS(patchPointLabels.size(), vector::zero);
//            scalarField ownSumW(patchPointLabels.size(), 0);
//            scalarField ownSumSqrW(patchPointLabels.size(), 0);

//            const labelList& nonGlobalPatchPoints =
//                procPatch.nonGlobalPatchPoints();

//            forAll(nonGlobalPatchPoints, pointI)
//            {
//                label curPatchPoint =
//                    nonGlobalPatchPoints[pointI];

//                label curPoint =
//                    patchPointLabels[curPatchPoint];

//                const labelList& curPointFaces = pointFaces[curPoint];

//                scalarField w(curPointFaces.size(), 0.0);

//                forAll (curPointFaces, faceI)
//                {
//                    label curFace = curPointFaces[faceI];

//                    w[faceI] = 1.0/mag
//                    (
//                        controlPoints()[curFace]
//                      - points[curPoint]
//                    );

//                    ownM[curPatchPoint] +=
//                        sqr(w[faceI])*sqr(controlPoints()[curFace]);

//                    ownS[curPatchPoint] +=
//                        sqr(w[faceI])*controlPoints()[curFace];
//                }

////#               include "emptyProcessorFaPatchPoints.H"

//                ownSumW[curPatchPoint] = sum(w);
//                ownSumSqrW[curPatchPoint] = sum(sqr(w));
//            }

//            // Parallel data exchange
//            {
//                OPstream toNeighbProc
//                (
//                    Pstream::blocking,
//                    procPatch.neighbProcNo(),
//                    ownM.size()*sizeof(symmTensor)
//                  + ownS.size()*sizeof(vector)
//                  + 4*ownSumW.size()*sizeof(scalar)
//                );

//                toNeighbProc << ownM << ownS << ownSumW << ownSumSqrW;
//            }

//            Field<symmTensor> ngbM
//            (
//                patchPointLabels.size(),
//                symmTensor::zero
//            );
//            vectorField ngbS
//            (
//                patchPointLabels.size(),
//                vector::zero
//            );
//            scalarField ngbSumW
//            (
//                patchPointLabels.size(),
//                0
//            );
//            scalarField ngbSumSqrW
//            (
//                patchPointLabels.size(),
//                0
//            );

//            {
//                IPstream fromNeighbProc
//                (
//                    Pstream::blocking,
//                    procPatch.neighbProcNo(),
//                    ngbM.size()*sizeof(symmTensor)
//                  + ngbS.size()*sizeof(vector)
//                  + 4*ngbSumW.size()*sizeof(scalar)
//                );

//                fromNeighbProc >> ngbM >> ngbS >> ngbSumW >> ngbSumSqrW;
//            }


//            forAll(nonGlobalPatchPoints, pointI)
//            {
//                label curPatchPoint =
//                    nonGlobalPatchPoints[pointI];

//                label curPoint =
//                    patchPointLabels[curPatchPoint];

//                label curNgbPoint = procPatch.neighbPoints()[curPatchPoint];
//                ownM[curPatchPoint] += ngbM[curNgbPoint];
//                ownS[curPatchPoint] += ngbS[curNgbPoint];
//                ownSumW[curPatchPoint] += ngbSumW[curNgbPoint];
//                ownSumSqrW[curPatchPoint] += ngbSumSqrW[curNgbPoint];

//                ownM[curPatchPoint] /= sqr(ownSumW[curPatchPoint]);
//                ownS[curPatchPoint] /= sqr(ownSumW[curPatchPoint]);

//                vector N = inv(ownM[curPatchPoint])&ownS[curPatchPoint];
//                N /= mag(N);

//                scalar p = (ownS[curPatchPoint]&N)
//                    /(ownSumSqrW[curPatchPoint]/sqr(ownSumW[curPatchPoint]));

//                displacement[curPoint] =
//                    pointsDisplacementDir()[curPoint]*
//                    (p - (points[curPoint]&N))
//                   /(pointsDisplacementDir()[curPoint]&N);
//            }
//        }
//    }


//    // Calculate displacement of global processor patch points
//    if (aMesh().globalData().nGlobalPoints() > 0)
//    {
//        const labelList& spLabels =
//            aMesh().globalData().sharedPointLabels();

//        Field<symmTensor> M(spLabels.size(), symmTensor::zero);
//        vectorField S(spLabels.size(), vector::zero);
//        scalarField sumW(spLabels.size(), 0);
//        scalarField sumSqrW(spLabels.size(), 0);

//        forAll (spLabels, pointI)
//        {
//            label curPoint = spLabels[pointI];

//            const labelList& curPointFaces = pointFaces[curPoint];

//            scalarField w(curPointFaces.size(), 0.0);

//            forAll (curPointFaces, faceI)
//            {
//                label curFace = curPointFaces[faceI];

//                w[faceI] = 1.0/mag
//                    (
//                        controlPoints()[curFace]
//                      - points[curPoint]
//                    );

//                M[pointI] =
//                    M[pointI]
//                  + sqr(w[faceI])*sqr(controlPoints()[curFace]);

//                    S[pointI] += sqr(w[faceI])*controlPoints()[curFace];
//            }

//            sumW[pointI] = sum(w);
//            sumSqrW[pointI] = sum(sqr(w));
//        }

//        const labelList& addr = aMesh().globalData().sharedPointAddr();
//        label nGlobalPoints = aMesh().globalData().nGlobalPoints();

//        symmTensorField gpM(nGlobalPoints, symmTensor::zero);
//        vectorField gpS(nGlobalPoints, vector::zero);
//        scalarField gpSumW(nGlobalPoints, 0.0);
//        scalarField gpSumSqrW(nGlobalPoints, 0.0);

//        forAll (addr, i)
//        {
//            gpM[addr[i]] += M[i];
//            gpS[addr[i]] += S[i];
//            gpSumW[addr[i]] += sumW[i];
//            gpSumSqrW[addr[i]] += sumSqrW[i];
//        }

//        combineReduce(gpM, plusEqOp<symmTensorField >());
//        combineReduce(gpS, plusEqOp<vectorField >());
//        combineReduce(gpSumW, plusEqOp<scalarField >());
//        combineReduce(gpSumSqrW, plusEqOp<scalarField >());

//        // Extract local data
//        forAll (addr, i)
//        {
//            M[i] = gpM[addr[i]];
//            S[i] = gpS[addr[i]];
//            sumW[i] = gpSumW[addr[i]];
//            sumSqrW[i] = gpSumSqrW[addr[i]];
//        }

//        M /= sqr(sumW);
//        S /= sqr(sumW);

//        vectorField N = inv(M)&S;
//        N /= mag(N);

//        scalarField p = (S&N)/(sumSqrW/sqr(sumW));

//        forAll(spLabels, pointI)
//        {
//            label curPoint = spLabels[pointI];

//            displacement[curPoint] =
//                pointsDisplacementDir()[curPoint]*
//                (p[pointI] - (points[curPoint]&N[pointI]))
//               /(pointsDisplacementDir()[curPoint]&N[pointI]);
//        }
//    }

//    return tdisplacement;
//}



tmp<vectorField> movingInterfacePatches::pointDisplacement(const scalarField& deltaH)
{
    Info<< "deltaHf:"
        << " sum local = " << gSum(mag(deltaH)) 
        << ", global = " << gSum(deltaH)
        << endl;

    const pointField& points = aMesh().patch().localPoints();
    const labelListList& pointFaces = aMesh().patch().pointFaces();

    controlPoints() += facesDisplacementDir()*deltaH;

        if (Pstream::master())
        {

            label patchID = aMesh().boundary().findPatchID("centerline");

            if (patchID != -1)
            {
                const labelList& eFaces =
                    aMesh().boundary()[patchID].edgeFaces();
                
    //            const labelList& eFaces =
    //                aMesh().boundary()[fixedPatchID].edgeFaces();

                const labelListList& fFaces = aMesh().patch().faceFaces();
                const vectorField& fCentres =
                    aMesh().areaCentres().internalField();

                forAll(eFaces, edgeI)
                {
                    const label& curFace = eFaces[edgeI];
                    const labelList& curFaceFaces = fFaces[curFace];

                    scalar H = 0.0;
                    label counter = 0;

                    forAll(curFaceFaces, faceI)
                    {
                        label index = findIndex(eFaces, curFaceFaces[faceI]);

                        if (index == -1)
                        {
                            H +=
                                facesDisplacementDir()[curFaceFaces[faceI]]
                            & (
                                    controlPoints()[curFaceFaces[faceI]]
                                - fCentres[curFaceFaces[faceI]]
                                );

                            counter++;
                        }
                    }

                    H /= counter;

                    controlPoints()[curFace] =
                        fCentres[curFace]
                    + facesDisplacementDir()[curFace]*H;
                }
            }
        }

    tmp<vectorField> tdisplacement
    (
        new vectorField
        (
            points.size(),
            vector::zero
        )
    );

    vectorField& displacement = tdisplacement();


    // Calculate displacement of internal points
    const vectorField& pointNormals = aMesh().pointAreaNormals();
    const edgeList& edges = aMesh().patch().edges();
    labelList internalPoints = aMesh().internalPoints();

    forAll (internalPoints, pointI)
    {
        label curPoint = internalPoints[pointI];

        const labelList& curPointFaces = pointFaces[curPoint];

        vectorField lsPoints(curPointFaces.size(), vector::zero);

        for (label i=0; i<curPointFaces.size(); i++)
        {
            label curFace = curPointFaces[i];

            lsPoints[i] = controlPoints()[curFace];
        }

        vectorField pointAndNormal =
            lsPlanePointAndNormal
            (
                lsPoints,
                points[curPoint],
                pointNormals[curPoint]
            );

        vector& P = pointAndNormal[0];
        vector& N = pointAndNormal[1];

        displacement[curPoint] =
            pointsDisplacementDir()[curPoint]
           *((P - points[curPoint])&N)
           /(pointsDisplacementDir()[curPoint]&N);
    }


    // Mirror control points
    FieldField<Field, vector> patchMirrorPoints(aMesh().boundary().size());

    forAll(patchMirrorPoints, patchI)
    {
        patchMirrorPoints.set
        (
            patchI,
            new vectorField
            (
                aMesh().boundary()[patchI].faPatch::size(),
                vector::zero
            )
        );

        vectorField N =
            aMesh().boundary()[patchI].ngbPolyPatchFaceNormals();

        const labelList peFaces =
            labelList::subList
            (
                aMesh().edgeOwner(),
                aMesh().boundary()[patchI].faPatch::size(),
                aMesh().boundary()[patchI].start()
            );

        const labelList& pEdges = aMesh().boundary()[patchI];

        vectorField peCentres(pEdges.size(), vector::zero);
        forAll(peCentres, edgeI)
        {
            peCentres[edgeI] =
                edges[pEdges[edgeI]].centre(points);
        }

        vectorField delta =
            vectorField(controlPoints(), peFaces)
          - peCentres;

        patchMirrorPoints[patchI] =
            peCentres + ((I - 2*N*N)&delta);
    }


    // Calculate displacement of boundary points
    labelList boundaryPoints = aMesh().boundaryPoints();

    const labelListList& edgeFaces = aMesh().patch().edgeFaces();
    const labelListList& pointEdges = aMesh().patch().pointEdges();

    forAll (boundaryPoints, pointI)
    {
        label curPoint = boundaryPoints[pointI];

        if (motionPointsMask()[curPoint] == 1)
        {
            // Calculating mirror points
            const labelList& curPointEdges = pointEdges[curPoint];

            vectorField mirrorPoints(2, vector::zero);

            label counter = -1;

            forAll (curPointEdges, edgeI)
            {
                label curEdge = curPointEdges[edgeI];

                if(edgeFaces[curEdge].size() == 1)
                {
                    label patchID = -1;
                    label edgeID = -1;
                    forAll(aMesh().boundary(), patchI)
                    {
                        const labelList& pEdges =
                            aMesh().boundary()[patchI];
                        label index = findIndex(pEdges, curEdge);
                        if (index != -1)
                        {
                            patchID = patchI;
                            edgeID = index;
                            break;
                        }
                    }

                    mirrorPoints[++counter] =
                        patchMirrorPoints[patchID][edgeID];
                }
            }

            // Calculating LS plane fit
            const labelList& curPointFaces = pointFaces[curPoint];

            vectorField lsPoints
            (
                curPointFaces.size() + mirrorPoints.size(),
                vector::zero
            );

            counter = -1;

            for (label i=0; i<curPointFaces.size(); i++)
            {
                label curFace = curPointFaces[i];

                lsPoints[++counter] = controlPoints()[curFace];
            }

            for (label i=0; i<mirrorPoints.size(); i++)
            {
                lsPoints[++counter] = mirrorPoints[i];
            }

            vectorField pointAndNormal =
                lsPlanePointAndNormal
                (
                    lsPoints,
                    points[curPoint],
                    pointNormals[curPoint]
                );

            vector& P = pointAndNormal[0];
            vector& N = pointAndNormal[1];

            displacement[curPoint] =
                pointsDisplacementDir()[curPoint]
               *((P - points[curPoint])&N)
               /(pointsDisplacementDir()[curPoint]&N);
        }
    }


    // Calculate displacement of axis point
    forAll (aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == wedgeFaPatch::typeName
        )
        {
            const wedgeFaPatch& wedgePatch =
                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

            if(wedgePatch.axisPoint() > -1)
            {
                label axisPoint = wedgePatch.axisPoint();

                displacement[axisPoint] =
                    pointsDisplacementDir()[axisPoint]
                   *(
                        pointsDisplacementDir()[axisPoint]
                       &(
                            controlPoints()[pointFaces[axisPoint][0]]
                          - points[axisPoint]
                        )
                    );
            }
        }
    }


    // Calculate displacement of processor patch points
    forAll (aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == processorFaPatch::typeName
        )
        {
            const processorFaPatch& procPatch =
                refCast<const processorFaPatch>(aMesh().boundary()[patchI]);

            const labelList& patchPointLabels =
                procPatch.pointLabels();

            FieldField<Field, vector> lsPoints(patchPointLabels.size());
            forAll(lsPoints, pointI)
            {
                lsPoints.set(pointI, new vectorField(0, vector::zero));
            }

            const labelList& nonGlobalPatchPoints =
                procPatch.nonGlobalPatchPoints();

            forAll(nonGlobalPatchPoints, pointI)
            {
                label curPatchPoint =
                    nonGlobalPatchPoints[pointI];

                label curPoint =
                    patchPointLabels[curPatchPoint];

                const labelList& curPointFaces = pointFaces[curPoint];

                lsPoints[curPatchPoint].setSize(curPointFaces.size());

                forAll(curPointFaces, faceI)
                {
                    label curFace = curPointFaces[faceI];

                    lsPoints[curPatchPoint][faceI] = controlPoints()[curFace];
                }

#               include "boundaryProcessorFaPatchPoints.H"
            }

            scalar lsPointsSize = 0;
            forAll(lsPoints, pointI)
            {
                lsPointsSize +=
                    2*lsPoints[pointI].size()*sizeof(vector);
            }

            // Parallel data exchange
            {
                OPstream toNeighbProc
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo(),
                    lsPointsSize
                );

                toNeighbProc << lsPoints;
            }

            FieldField<Field, vector> ngbLsPoints(patchPointLabels.size());

            {
                IPstream fromNeighbProc
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo(),
                    lsPointsSize
                );

                fromNeighbProc >> ngbLsPoints;
            }

            forAll(nonGlobalPatchPoints, pointI)
            {
                label curPatchPoint =
                    nonGlobalPatchPoints[pointI];

                label curPoint =
                    patchPointLabels[curPatchPoint];

                label curNgbPoint = procPatch.neighbPoints()[curPatchPoint];

                vectorField allLsPoints
                (
                    lsPoints[curPatchPoint].size()
                  + ngbLsPoints[curNgbPoint].size(),
                    vector::zero
                );

                label counter = -1;
                forAll(lsPoints[curPatchPoint], pointI)
                {
                    allLsPoints[++counter] = lsPoints[curPatchPoint][pointI];
                }
                forAll(ngbLsPoints[curNgbPoint], pointI)
                {
                    allLsPoints[++counter] = ngbLsPoints[curNgbPoint][pointI];
                }

                vectorField pointAndNormal =
                    lsPlanePointAndNormal
                    (
                        allLsPoints,
                        points[curPoint],
                        pointNormals[curPoint]
                    );

                vector& P = pointAndNormal[0];
                vector& N = pointAndNormal[1];

                if (motionPointsMask()[curPoint] != 0)
                {
                    displacement[curPoint] =
                        pointsDisplacementDir()[curPoint]
                       *((P - points[curPoint])&N)
                       /(pointsDisplacementDir()[curPoint]&N);
                }
            }
        }
    }


    // Calculate displacement of global processor patch points
    if (aMesh().globalData().nGlobalPoints() > 0)
    {
        const labelList& spLabels =
            aMesh().globalData().sharedPointLabels();

        const labelList& addr = aMesh().globalData().sharedPointAddr();

        for (label k=0; k<aMesh().globalData().nGlobalPoints(); k++)
        {
            List<List<vector> > procLsPoints(Pstream::nProcs());

            label curSharedPointIndex = findIndex(addr, k);

            if (curSharedPointIndex != -1)
            {
                label curPoint = spLabels[curSharedPointIndex];

                const labelList& curPointFaces = pointFaces[curPoint];

                procLsPoints[Pstream::myProcNo()] =
                    List<vector>(curPointFaces.size());

                forAll (curPointFaces, faceI)
                {
                    label curFace = curPointFaces[faceI];

                    procLsPoints[Pstream::myProcNo()][faceI] =
                        controlPoints()[curFace];
                }
            }

            Pstream::gatherList(procLsPoints);
            Pstream::scatterList(procLsPoints);

            if (curSharedPointIndex != -1)
            {
                label curPoint = spLabels[curSharedPointIndex];

                label nAllPoints = 0;
                forAll(procLsPoints, procI)
                {
                    nAllPoints += procLsPoints[procI].size();
                }

                vectorField allPoints(nAllPoints, vector::zero);

                label counter = 0;
                forAll(procLsPoints, procI)
                {
                    forAll(procLsPoints[procI], pointI)
                    {
                        allPoints[counter++] =
                            procLsPoints[procI][pointI];
                    }
                }

                vectorField pointAndNormal =
                    lsPlanePointAndNormal
                    (
                        allPoints,
                        points[curPoint],
                        pointNormals[curPoint]
                    );

                const vector& P = pointAndNormal[0];
                const vector& N = pointAndNormal[1];

                displacement[curPoint] =
                    pointsDisplacementDir()[curPoint]
                   *((P - points[curPoint])&N)
                   /(pointsDisplacementDir()[curPoint]&N);
            }
        }
    }

    Info<< "displacement:"
        << " sum local = " << gSum(mag(displacement)) 
        << ", global = " << gSum(displacement)
        << endl;

    return tdisplacement;
}




//tmp<vectorField> movingInterfacePatches::pointDisplacement(const scalarField& deltaH)
//{
//    const pointField& points = aMesh().patch().localPoints();
//    const labelListList& pointFaces = aMesh().patch().pointFaces();

//    controlPoints() += facesDisplacementDir()*deltaH;

//        if (Pstream::master())
//        {

//            label patchID = aMesh().boundary().findPatchID("centerline");

//            const labelList& eFaces =
//                aMesh().boundary()[patchID].edgeFaces();

////            const labelList& eFaces =
////                aMesh().boundary()[fixedPatchID].edgeFaces();

//            const labelListList& fFaces = aMesh().patch().faceFaces();
//            const vectorField& fCentres =
//                aMesh().areaCentres().internalField();

//            forAll(eFaces, edgeI)
//            {
//                const label& curFace = eFaces[edgeI];
//                const labelList& curFaceFaces = fFaces[curFace];

//                scalar H = 0.0;
//                label counter = 0;

//                forAll(curFaceFaces, faceI)
//                {
//                    label index = findIndex(eFaces, curFaceFaces[faceI]);

//                    if (index == -1)
//                    {
//                        H +=
//                            facesDisplacementDir()[curFaceFaces[faceI]]
//                          & (
//                                controlPoints()[curFaceFaces[faceI]]
//                              - fCentres[curFaceFaces[faceI]]
//                            );

//                        counter++;
//                    }
//                }

//                H /= counter;

//                controlPoints()[curFace] =
//                    fCentres[curFace]
//                  + facesDisplacementDir()[curFace]*H;
//            }
//        }

//    tmp<vectorField> tdisplacement
//    (
//        new vectorField
//        (
//            points.size(),
//            vector::zero
//        )
//    );

//    vectorField& displacement = tdisplacement();


//    // Calculate displacement of internal points
//    const vectorField& pointNormals = aMesh().pointAreaNormals();
//    const edgeList& edges = aMesh().patch().edges();
//    labelList internalPoints = aMesh().internalPoints();

//    forAll (internalPoints, pointI)
//    {
//        label curPoint = internalPoints[pointI];

//        const labelList& curPointFaces = pointFaces[curPoint];

//        vectorField lsPoints(curPointFaces.size(), vector::zero);

//        for (label i=0; i<curPointFaces.size(); i++)
//        {
//            label curFace = curPointFaces[i];

//            lsPoints[i] = controlPoints()[curFace];
//        }

//        vectorField pointAndNormal =
//            lsPlanePointAndNormal
//            (
//                lsPoints,
//                points[curPoint],
//                pointNormals[curPoint]
//            );

//        vector& P = pointAndNormal[0];
//        vector& N = pointAndNormal[1];

//        displacement[curPoint] =
//            pointsDisplacementDir()[curPoint]
//           *((P - points[curPoint])&N)
//           /(pointsDisplacementDir()[curPoint]&N);
//    }


//    // Mirror control points
//    FieldField<Field, vector> patchMirrorPoints(aMesh().boundary().size());

//    forAll(patchMirrorPoints, patchI)
//    {
//        patchMirrorPoints.set
//        (
//            patchI,
//            new vectorField
//            (
//                aMesh().boundary()[patchI].faPatch::size(),
//                vector::zero
//            )
//        );

//        vectorField N =
//            aMesh().boundary()[patchI].ngbPolyPatchFaceNormals();

//        const labelList peFaces =
//            labelList::subList
//            (
//                aMesh().edgeOwner(),
//                aMesh().boundary()[patchI].faPatch::size(),
//                aMesh().boundary()[patchI].start()
//            );

//        const labelList& pEdges = aMesh().boundary()[patchI];

//        vectorField peCentres(pEdges.size(), vector::zero);
//        forAll(peCentres, edgeI)
//        {
//            peCentres[edgeI] =
//                edges[pEdges[edgeI]].centre(points);
//        }

//        vectorField delta =
//            vectorField(controlPoints(), peFaces)
//          - peCentres;

//        patchMirrorPoints[patchI] =
//            peCentres + ((I - 2*N*N)&delta);
//    }


//    // Calculate displacement of boundary points
//    labelList boundaryPoints = aMesh().boundaryPoints();

//    const labelListList& edgeFaces = aMesh().patch().edgeFaces();
//    const labelListList& pointEdges = aMesh().patch().pointEdges();

//    forAll (boundaryPoints, pointI)
//    {
//        label curPoint = boundaryPoints[pointI];

//        if (motionPointsMask()[curPoint] == 1)
//        {
//            // Calculating mirror points
//            const labelList& curPointEdges = pointEdges[curPoint];

//            vectorField mirrorPoints(2, vector::zero);

//            label counter = -1;

//            forAll (curPointEdges, edgeI)
//            {
//                label curEdge = curPointEdges[edgeI];

//                if(edgeFaces[curEdge].size() == 1)
//                {
//                    label patchID = -1;
//                    label edgeID = -1;
//                    forAll(aMesh().boundary(), patchI)
//                    {
//                        const labelList& pEdges =
//                            aMesh().boundary()[patchI];
//                        label index = findIndex(pEdges, curEdge);
//                        if (index != -1)
//                        {
//                            patchID = patchI;
//                            edgeID = index;
//                            break;
//                        }
//                    }

//                    mirrorPoints[++counter] =
//                        patchMirrorPoints[patchID][edgeID];
//                }
//            }

//            // Calculating LS plane fit
//            const labelList& curPointFaces = pointFaces[curPoint];

//            vectorField lsPoints
//            (
//                curPointFaces.size() + mirrorPoints.size(),
//                vector::zero
//            );

//            counter = -1;

//            for (label i=0; i<curPointFaces.size(); i++)
//            {
//                label curFace = curPointFaces[i];

//                lsPoints[++counter] = controlPoints()[curFace];
//            }

//            for (label i=0; i<mirrorPoints.size(); i++)
//            {
//                lsPoints[++counter] = mirrorPoints[i];
//            }

//            vectorField pointAndNormal =
//                lsPlanePointAndNormal
//                (
//                    lsPoints,
//                    points[curPoint],
//                    pointNormals[curPoint]
//                );

//            vector& P = pointAndNormal[0];
//            vector& N = pointAndNormal[1];

//            displacement[curPoint] =
//                pointsDisplacementDir()[curPoint]
//               *((P - points[curPoint])&N)
//               /(pointsDisplacementDir()[curPoint]&N);
//        }
//    }


////    // Calculate displacement of axis point
////    forAll (aMesh().boundary(), patchI)
////    {
////        if
////        (
////            aMesh().boundary()[patchI].type()
////         == wedgeFaPatch::typeName
////        )
////        {
////            const wedgeFaPatch& wedgePatch =
////                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

////            if(wedgePatch.axisPoint() > -1)
////            {
////                label axisPoint = wedgePatch.axisPoint();

////                displacement[axisPoint] =
////                    pointsDisplacementDir()[axisPoint]
////                   *(
////                        pointsDisplacementDir()[axisPoint]
////                       &(
////                            controlPoints()[pointFaces[axisPoint][0]]
////                          - points[axisPoint]
////                        )
////                    );
////            }
////        }
////    }


//    // Calculate displacement of processor patch points
//    forAll (aMesh().boundary(), patchI)
//    {
//        if
//        (
//            aMesh().boundary()[patchI].type()
//         == processorFaPatch::typeName
//        )
//        {
//            const processorFaPatch& procPatch =
//                refCast<const processorFaPatch>(aMesh().boundary()[patchI]);

//            const labelList& patchPointLabels =
//                procPatch.pointLabels();

//            FieldField<Field, vector> lsPoints(patchPointLabels.size());
//            forAll(lsPoints, pointI)
//            {
//                lsPoints.set(pointI, new vectorField(0, vector::zero));
//            }

//            const labelList& nonGlobalPatchPoints =
//                procPatch.nonGlobalPatchPoints();

//            forAll(nonGlobalPatchPoints, pointI)
//            {
//                label curPatchPoint =
//                    nonGlobalPatchPoints[pointI];

//                label curPoint =
//                    patchPointLabels[curPatchPoint];

//                const labelList& curPointFaces = pointFaces[curPoint];

//                lsPoints[curPatchPoint].setSize(curPointFaces.size());

//                forAll(curPointFaces, faceI)
//                {
//                    label curFace = curPointFaces[faceI];

//                    lsPoints[curPatchPoint][faceI] = controlPoints()[curFace];
//                }

//#               include "boundaryProcessorFaPatchPoints.H"
//            }

//            scalar lsPointsSize = 0;
//            forAll(lsPoints, pointI)
//            {
//                lsPointsSize +=
//                    2*lsPoints[pointI].size()*sizeof(vector);
//            }

//            // Parallel data exchange
//            {
//                OPstream toNeighbProc
//                (
//                    Pstream::nonBlocking,
//                    procPatch.neighbProcNo(),
//                    lsPointsSize
//                );

//                toNeighbProc << lsPoints;
//            }

//            FieldField<Field, vector> ngbLsPoints(patchPointLabels.size());

//            {
//                IPstream fromNeighbProc
//                (
//                    Pstream::blocking,
//                    procPatch.neighbProcNo(),
//                    lsPointsSize*3
//                );

//                fromNeighbProc >> ngbLsPoints;
//            }

//            forAll(nonGlobalPatchPoints, pointI)
//            {
//                label curPatchPoint =
//                    nonGlobalPatchPoints[pointI];

//                label curPoint =
//                    patchPointLabels[curPatchPoint];

//                label curNgbPoint = procPatch.neighbPoints()[curPatchPoint];

//                vectorField allLsPoints
//                (
//                    lsPoints[curPatchPoint].size()
//                  + ngbLsPoints[curNgbPoint].size(),
//                    vector::zero
//                );

//                label counter = -1;
//                forAll(lsPoints[curPatchPoint], pointI)
//                {
//                    allLsPoints[++counter] = lsPoints[curPatchPoint][pointI];
//                }
//                forAll(ngbLsPoints[curNgbPoint], pointI)
//                {
//                    allLsPoints[++counter] = ngbLsPoints[curNgbPoint][pointI];
//                }

//                vectorField pointAndNormal =
//                    lsPlanePointAndNormal
//                    (
//                        allLsPoints,
//                        points[curPoint],
//                        pointNormals[curPoint]
//                    );

//                vector& P = pointAndNormal[0];
//                vector& N = pointAndNormal[1];

//                if (motionPointsMask()[curPoint] != 0)
//                {
//                    displacement[curPoint] =
//                        pointsDisplacementDir()[curPoint]
//                       *((P - points[curPoint])&N)
//                       /(pointsDisplacementDir()[curPoint]&N);
//                }
//            }
//        }
//    }


//    // Calculate displacement of global processor patch points
//    if (aMesh().globalData().nGlobalPoints() > 0)
//    {
//        const labelList& spLabels =
//            aMesh().globalData().sharedPointLabels();

//        const labelList& addr = aMesh().globalData().sharedPointAddr();

//        for (label k=0; k<aMesh().globalData().nGlobalPoints(); k++)
//        {
//            List<List<vector> > procLsPoints(Pstream::nProcs());

//            label curSharedPointIndex = findIndex(addr, k);

//            if (curSharedPointIndex != -1)
//            {
//                label curPoint = spLabels[curSharedPointIndex];

//                const labelList& curPointFaces = pointFaces[curPoint];

//                procLsPoints[Pstream::myProcNo()] =
//                    List<vector>(curPointFaces.size());

//                forAll (curPointFaces, faceI)
//                {
//                    label curFace = curPointFaces[faceI];

//                    procLsPoints[Pstream::myProcNo()][faceI] =
//                        controlPoints()[curFace];
//                }
//            }

//            Pstream::gatherList(procLsPoints);
//            Pstream::scatterList(procLsPoints);

//            if (curSharedPointIndex != -1)
//            {
//                label curPoint = spLabels[curSharedPointIndex];

//                label nAllPoints = 0;
//                forAll(procLsPoints, procI)
//                {
//                    nAllPoints += procLsPoints[procI].size();
//                }

//                vectorField allPoints(nAllPoints, vector::zero);

//                label counter = 0;
//                forAll(procLsPoints, procI)
//                {
//                    forAll(procLsPoints[procI], pointI)
//                    {
//                        allPoints[counter++] =
//                            procLsPoints[procI][pointI];
//                    }
//                }

//                vectorField pointAndNormal =
//                    lsPlanePointAndNormal
//                    (
//                        allPoints,
//                        points[curPoint],
//                        pointNormals[curPoint]
//                    );

//                const vector& P = pointAndNormal[0];
//                const vector& N = pointAndNormal[1];

//                displacement[curPoint] =
//                    pointsDisplacementDir()[curPoint]
//                   *((P - points[curPoint])&N)
//                   /(pointsDisplacementDir()[curPoint]&N);
//            }
//        }
//    }

//    return tdisplacement;
//}

Foam::tmp<vectorField>
Foam::movingInterfacePatches::lsPlanePointAndNormal
(
    const vectorField& points,
    const vector& origin,
    const vector& axis
) const
{
    // LS in local CS
    vector dir = (points[0] - origin);
    dir -= axis*(axis&dir);
    dir /= mag(dir);
    coordinateSystem cs("cs", origin, axis, dir);

    vectorField localPoints = cs.localPosition(points);
    scalarField W = 1.0/(mag(points - origin) + SMALL);

    scalarRectangularMatrix M
    (
        points.size(),
        3,
        0.0
    );

    for (label i=0; i<localPoints.size(); i++)
    {
        M[i][0] = localPoints[i].x();
        M[i][1] = localPoints[i].y();
        M[i][2] = 1.0;
    }

    scalarSquareMatrix MtM(3, 0.0);
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

    scalarField MtR(3, 0);
    for (label i = 0; i < MtR.size(); i++)
    {
        for (label j = 0; j < M.n(); j++)
        {
            MtR[i] += M[j][i]*localPoints[j].z()*W[j];
        }
    }

    scalarSquareMatrix::LUsolve(MtM, MtR);

    vector n0 = vector(-MtR[0], -MtR[1], 1);
    n0 = cs.globalVector(n0);
    n0 /= mag(n0);

    vector p0 = vector(0, 0, MtR[2]);
    p0 = cs.globalPosition(p0);

    tmp<vectorField> pointAndNormal
    (
        new vectorField(2, vector::zero)
    );

    pointAndNormal()[0] = p0;
    pointAndNormal()[1] = n0;

    return pointAndNormal;
}


// ************************************************************************* //
