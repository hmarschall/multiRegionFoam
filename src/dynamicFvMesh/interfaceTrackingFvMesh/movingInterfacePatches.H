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

Class
    movingInterfacePatches

Description
    Move surface patches for ALE interface-tracking

    TODO:
    Trigger mesh motion if threshold for surface patch motion is hit. These 
    should then be collected over all surfaces. If one of the surfaces triggers, 
    do mesh motion.

Author
    Holger Marschall <holger.marschall@tu-darmstadt.de> 

SourceFiles
    movingInterfacePatches.C

\*---------------------------------------------------------------------------*/

#ifndef movingInterfacePatches_H
#define movingInterfacePatches_H

#include "Switch.H"
#include "polyMesh.H"
#include "dynamicFvMesh.H"
#include "demandDrivenData.H"
#include "globalPolyPatch.H"
#include "interfaceToInterfaceMapping.H"

#include "faCFD.H"
#include "fvCFD.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyTopoChanger;

/*---------------------------------------------------------------------------*\
                        Class movingInterfacePatches Declaration
\*---------------------------------------------------------------------------*/

class movingInterfacePatches
{
    // Private data

        //- Name
        const word name_;

        //- Reference to the mesh
        const dynamicFvMesh& mesh_;

        //- Finite area mesh
        autoPtr<faMesh> aMeshPtr_;

        //- Motion sub-dictionary
        const dictionary& motionDict_;

        //- Markup field for points.  Moving points marked with 1
//        mutable scalarField* movingPointsMaskPtr_;

        //- Surface patch index
        label surfacePatchID_;

        //- Reference to the patch
        const fvPatch& surfacePatch_;

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

        //- Interface smoothing
        Switch smoothing_;

        //- Smoothing frequency
        label smoothingFrequency_;

        //- Smoothing tolerance
        scalar smoothingTolerance_;

        //- Rigid free-surface
        Switch rigidFreeSurface_;

        //- Interface? (patch and shadow patch to move)
        bool isInterface_;

        //- Move surface and mesh?
        bool motionTrigger_;

        //- Limit above which interface is deformed
        dimensionedScalar interfaceDeformationLimit_;

        //- Current interface tracking time index
        mutable label timeIndex_;

        //- Swept Volume Correction
        scalarField sweptVolCorrOld_;

        //- Reset flux frequency
        label resetFluxFrequency_;

        //- Neighbour region (aka mesh) name
        word neighbourRegionName_;

        //- Neighbour (aka shadow) patch name
        word neighbourPatchName_;

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

        //- Global poly patch
        //  This is a copy of a patch which is reconstructed in full on each
        //  processor. This is convenient for interpolating between two patches
        //  in parallel
        //- A-side
        mutable autoPtr<globalPolyPatch> globalPatchPtr_;
        //- B-side
        mutable autoPtr<globalPolyPatch> globalNbrPatchPtr_;

        //- Mapping object for passing fields between the pairs of interfaces
        //  e.g. this could be GGI or a direct map or other
        mutable autoPtr<interfaceToInterfaceMapping> interfaceToInterfacePtr_;

        //- Interpolator update frequency
        label interpolatorUpdateFrequency_;


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

            //- Make global poly patch
            void makeGlobalPatches() const;

            //- Clear global poly patch
            void clearGlobalPatches() const;

            //- Make interface-to-interface
            void makeInterfaceToInterface() const;

        //- Disallow default bitwise copy construct
        movingInterfacePatches(const movingInterfacePatches&);

        //- Disallow default bitwise assignment
        void operator=(const movingInterfacePatches&);

        //- Calculate moving mask
//        void calcMovingMask() const;

        //- Return moving points mask
//        const scalarField& movingPointsMask() const;

        //- Clear moving points mask
//        void clearPointMask();


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
    TypeName("movingInterfacePatches");

    // Constructors

        //- Construct from  dictionary
        movingInterfacePatches
        (
            const word& name,
            const dynamicFvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~movingInterfacePatches();


    // Member Functions

        //- Return reference to finite volume mesh
        const dynamicFvMesh& mesh() const
        {
            return mesh_;
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

        const label& patchID() const
        {
            return surfacePatchID_;
        }

        const fvPatch& patch() const
        {
            return surfacePatch_;
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

        //- Correct the point normal vectors
        void correctPointNormals();

        //- Motion direction switch
        bool normalMotionDir() const
        {
            return normalMotionDir_;
        }

        //- Neighbour mesh
        const Foam::dynamicFvMesh& nbrMesh() const;

        //- Neighbour patch
        const Foam::fvPatch& nbrPatch() const;

        //- Return const access to the interface-to-interface object
        const interfaceToInterfaceMapping& interfaceToInterface() const;

        //- Global poly patch

        //-- Return const-access to the global poly patch A
        const globalPolyPatch& globalPatch() const;

        //-- Return const-access to the global poly patch B
        const globalPolyPatch& globalNbrPatch() const;

        //- Clear control points
        void clearControlPoints()
        {
            controlPointsPtr_.clear();
        }

        const bool& isInterface()
        {
            return isInterface_;
        }

        const bool& moving()
        {
            return motionTrigger_;
        }

        //- Update interpolator
        virtual void updateInterpolatorAndGlobalPatches();

        //- Write VTK surface mesh
        void writeVTK() const;

        //- Write VTK surface control points
        void writeVTKControlPoints();

        //- Return surface point motion
        virtual tmp<vectorField> surfacePointDisplacement();

        //- Return shadow point motion
        virtual tmp<vectorField> shadowPointDisplacement
        (
            const vectorField& displacement
        );

        //- Smooth surface mesh using Laplacian smoothing
        virtual tmp<vectorField> smoothSurfaceMesh();

        //- Update topology
        void updateTopology();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
