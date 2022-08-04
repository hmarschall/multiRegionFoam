/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description

\*---------------------------------------------------------------------------*/

#include "faMesh.H"
#include "foamTime.H"
#include "areaFields.H"
#include "edgeFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::faMesh::movePoints() const
// Foam::faMesh::movePoints(const vectorField& newPoints)
{
    // Grab point motion from polyMesh
    const vectorField& newPoints = mesh().allPoints();

    // Grab old time areas if the time has been incremented
    if (curTimeIndex_ < time().timeIndex())
    {
        if (S00Ptr_ && S0Ptr_)
        {
            Info<< "Copy old-old S" << endl;
            *S00Ptr_ = *S0Ptr_;
        }

        if (S0Ptr_)
        {
            Info<< "Copy old S" << endl;
            *S0Ptr_ = S();
        }
        else
        {
            if (debug)
            {
                InfoIn("bool faMesh::movePoints() const")
                    << "Creating old cell volumes." << endl;
            }

            S0Ptr_ = new DimensionedField<scalar, areaMesh>
            (
                IOobject
                (
                    "S0",
                    time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                S()
            );
        }

        curTimeIndex_ = time().timeIndex();
    }

    clearGeomNotAreas();
//    clearOut();

    pointAreaNormals();

    // To satisfy the motion interface for MeshObject, const cast is needed
    // HJ, 5/Aug/2011
//    if (patchPtr_)
//    {
//        patchPtr_->movePoints(newPoints);
//    }

    // Move boundary points
    const_cast<faBoundaryMesh&>(boundary_).movePoints(newPoints);

    // Move interpolation
    const edgeInterpolation& cei = *this;
    const_cast<edgeInterpolation&>(cei).edgeInterpolation::movePoints();

    clearGeomNotAreas();

//     // Create global mesh data
//     if (Pstream::parRun())
//     {
//         globalData();
//     }

    // Calculate topology for the patches (processor-processor comms etc.)
    const_cast<faBoundaryMesh&>(boundary_).updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    const_cast<faBoundaryMesh&>(boundary_).calcGeometry();

    // Fluxes were dummy?  HJ, 28/Jul/2011

    return true;



//     moving_ = true;

//     // Grab old time areas if the time has been incremented
//     if (curMotionTimeIndex_ < operator()().time().timeIndex())
//     {
//         if (S00Ptr_ && S0Ptr_)
//         {
//             *S00Ptr_ = *S0Ptr_;
//         }

//         if (S0Ptr_)
//         {
//             *S0Ptr_ = S();
//         }
//         else
//         {
//             S0Ptr_ = new DimensionedField<scalar, areaMesh>
//             (
//                 IOobject
//                 (
//                     "S0",
//                     time().timeName(),
//                     mesh_,
//                     IOobject::NO_READ,
//                     IOobject::NO_WRITE,
//                     false
//                 ),
//                 S()
//             );
//         }

//         curMotionTimeIndex_ = operator()().time().timeIndex();
//     }

//     clearGeomNotAreas();

//     pointAreaNormals();

//     patch().movePoints(newPoints);
//     boundary_.movePoints(newPoints);
//     edgeInterpolation::movePoints();

//     clearGeomNotAreas();

//     tmp<scalarField> tresult(new scalarField(nEdges(), 0.0));

//     return tresult;
}

// ************************************************************************* //
