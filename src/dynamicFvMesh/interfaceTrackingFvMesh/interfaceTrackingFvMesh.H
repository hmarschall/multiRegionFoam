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

Class
    Foam::interfaceTrackingFvMesh

Description
    The interfaceTrackingFvMesh

SourceFiles
    interfaceTrackingFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceTrackingFvMesh_H
#define interfaceTrackingFvMesh_H

//#include "dynamicMotionSolverFvMesh.H"
#include "topoChangerFvMesh.H"
//#include "dynamicFvMesh.H"
#include "motionSolver.H"
#include "regIOobject.H"
#include "faCFD.H"
#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "singlePhaseTransportModel.H"
#include "demandDrivenData.H"
//#include "correctPointsFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class interfaceTrackingFvMesh Declaration
\*---------------------------------------------------------------------------*/

class interfaceTrackingFvMesh
:
    public topoChangerFvMesh//dynamicMotionSolverFvMesh
{
    // Private Data

        // Dynamic motion solver
        autoPtr<motionSolver> motionPtr_;

        //- Finite area mesh
        autoPtr<faMesh> aMeshPtr_;

        //- Motion dictionary
        dictionary motionDict_;

        //- Surface patches which move
        wordList movingSurfacePatches_;

        //- Surface patch index
        label surfacePatchID_;

        //- Surface patches which do not move
        wordList fixedSurfacePatches_;

        //- Surface faPatch-es where wave should not reflect
        wordList nonReflectingFreeSurfacePatches_;

        //- Surface patches for which point normals must be corrected
        wordList pointNormalsCorrectionPatches_;

        //- Is the surface points displacement direction
        //  parallel with surface point normals?
        Switch normalMotionDir_;

        //- Surface points displacement direction
        //  if not normal motion direction
        vector motionDir_;

        //- Interface smoothing at the beginning of time step
        Switch smoothing_;

        //- Rigid free-surface
        Switch rigidFreeSurface_;

        //- Current interface tracking time index
        label timeIndex_;

        //- Swept Volume Correction
        scalarField sweptVolCorrOld_;

        //- Reset flux frequency
        label resetFluxFrequency_;

        //- Points which are attached to the free-surface A side faces
        //  and which defines the free-surface shape
        autoPtr<vectorIOField> controlPointsPtr_;

        //- Total displacement of free-surface points
        //  in one time step
        autoPtr<vectorIOField> totalDisplacementPtr_;

        //- Field which additionally determines
        //  the motion of free-surface points
//        mutable labelList* motionPointsMaskPtr_;
        autoPtr<scalarField> motionPointsMaskPtr_;

        //- Displacement direction of free-surface points
        autoPtr<vectorField> pointsDisplacementDirPtr_;

        //- Displacement direction of free-surface control points
        autoPtr<vectorField> facesDisplacementDirPtr_;

        //- Contact angle
//        mutable areaScalarField* contactAnglePtr_;

    // Private Member Functions

        // Make demand-driven data

            //- Create control points
            void makeControlPoints();

            //- Create total displacement field
            void makeTotalDisplacement();

            //- Create motion points mask
            void makeMotionPointsMask();

            //- Create points and control point motion direction
            void makeDirections();

            //- Create contact angle
//            void makeContactAngle();

        //- No copy construct
        interfaceTrackingFvMesh(const interfaceTrackingFvMesh&) = delete;

        //- No copy assignment
        void operator=(const interfaceTrackingFvMesh&) = delete;

        //- Initialize data
        void initializeData();

        //- Update control points end displacement directions
        void updateDisplacementDirections();

        //- Initialize control points position
        void initializeControlPointsPosition();

        //- Calculate free surface points displacement
        //  for given new control points position
        tmp<vectorField> pointDisplacement
        (
            const scalarField& deltaH
        );

        //- Calc least sqare plane point and normal
        tmp<vectorField> lsPlanePointAndNormal
        (
            const vectorField& points,
            const vector& origin,
            const vector& axis
        ) const;

        //- Correct free surface point displacement next to fixed contact line
        void correctPointDisplacement
        (
            const scalarField& sweptVolCorr,
            vectorField& displacement
        );

public:

    //- Runtime type information
    TypeName("interfaceTrackingFvMesh");


    // Constructors

        //- Construct from IOobject
        interfaceTrackingFvMesh(const IOobject& io);


    //- Destructor
    ~interfaceTrackingFvMesh();


    // Member Functions

        fvMesh& mesh()
        {
            return *this;
        }

        const fvMesh& mesh() const
        {
            return *this;
        }

        //- Return reference to finite area mesh
        faMesh& aMesh()
        {
            return aMeshPtr_();
        }

        //- Return reference to finite area mesh
        const faMesh& aMesh() const
        {
            return aMeshPtr_();
        }

        const label& surfacePatchID() const
        {
            return surfacePatchID_;
        }

        //- Rigid free-surface
        const Switch& rigidFreeSurface() const
        {
            return rigidFreeSurface_;
        }

        //- Rigid free-surface
        Switch& rigidFreeSurface()
        {
            return rigidFreeSurface_;
        }

//        //- Correct contact line point normals
//        const Switch& correctContactLineNormals() const
//        {
//            return correctContactLineNormals_;
//        }

//        //- Correct contact line point normals
//        Switch& correctContactLineNormals()
//        {
//            return correctContactLineNormals_;
//        }

        //- Return constant reference to velocity field
        const volVectorField& U() const;

        //- Return constant reference to pressure field
        const volScalarField& p() const;

        //- Return constant reference to velocity field
        const surfaceScalarField& phi() const;

        //- Return control points
        vectorField& controlPoints();

        //- Return reference to motion points mask field
//        labelList& motionPointsMask();
        scalarField& motionPointsMask();

        //- Return reference to point displacement direction field
        vectorField& pointsDisplacementDir();

        //- Return reference to control points displacement direction field
        vectorField& facesDisplacementDir();

        //- Return total points displacement
        vectorField& totalDisplacement();

        //- Motion direction swithc
        bool normalMotionDir() const
        {
            return normalMotionDir_;
        }

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Clear control points
        void clearControlPoints()
        {
            controlPointsPtr_.clear();
        }

        //- Write VTK surface mesh
        void writeVTK() const;

        //- Write VTK surface control points
        void writeVTKControlPoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
