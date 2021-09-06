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

#include "interfaceTrackingFvMesh.H"
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

Foam::tmp<Foam::vectorField>
Foam::interfaceTrackingFvMesh::pointDisplacement(const scalarField& deltaH)
{
    const pointField& points = aMesh().patch().localPoints();
    const labelListList& pointFaces = aMesh().patch().pointFaces();

    const labelList faceCells = 
        mesh().boundary()[surfacePatchID()].patch().faceCells();

    controlPoints() += facesDisplacementDir()*deltaH;

    // Correct curvature at axis
    if (Pstream::master())
    {
        label patchID = aMesh().boundary().findPatchID("centerline");

        if (patchID != -1)
        {
            const labelList& eFaces =
                aMesh().boundary()[patchID].edgeFaces();

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

    // Correct controPoints next to fixed patches
//    {
        forAll(fixedSurfacePatches_, patchI)
        {
            label fixedPatchID =
                aMesh().boundary().findPatchID
                (
                    fixedSurfacePatches_[patchI]
                );

            if(fixedPatchID == -1)
            {
                FatalErrorIn("freeSurface::freeSurface(...)")
                    << "Wrong faPatch name in the fixedSurfacePatches list"
                        << " defined in the freeSurfaceProperties dictionary"
                        << abort(FatalError);
            }

            const labelList& eFaces =
                aMesh().boundary()[fixedPatchID].edgeFaces();

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
//    }

    // Calculate displacement of internal points
    tmp<vectorField> tdisplacement
    (
        new vectorField
        (
            points.size(),
            vector::zero
        )
    );

    vectorField& displacement = tdisplacement();

    forAll (pointFaces, pointI)
    {
        scalar weightsSum = 0.0;
        const labelList& curPointFaces = pointFaces[pointI];

        forAll (curPointFaces, faceI)
        {
            label curFace = curPointFaces[faceI];

            scalar weight = 1.0/mag
            (
                points[pointI]
              - controlPoints()[curFace]
            );

            displacement[pointI] += weight*controlPoints()[curFace];

            weightsSum += weight;
        }

        displacement[pointI] /= weightsSum;

        displacement[pointI] -= points[pointI];
    }

    displacement = motionPointsMask()*
        (pointsDisplacementDir()&displacement)*
        pointsDisplacementDir();

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

    return tdisplacement;
}


Foam::tmp<vectorField>
Foam::interfaceTrackingFvMesh::lsPlanePointAndNormal
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
