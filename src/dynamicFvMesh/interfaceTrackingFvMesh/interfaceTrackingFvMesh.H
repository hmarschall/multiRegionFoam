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
    Dynamic mesh class for ALE interface-tracking supporting multiple 
    moving interface patches, where interface refers to calculated patch motion 
    and consistent neigbouring mesh patch motion (shadow motion)

    Contact line support still to be implemented.

Author
    Holger Marschall <holger.marschall@tu-darmstadt.de> 

SourceFiles
    interfaceTrackingFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceTrackingFvMesh_H
#define interfaceTrackingFvMesh_H

//#include "dynamicMotionSolverFvMesh.H"
//#include "dynamicFvMesh.H"
#include "topoChangerFvMesh.H"
#include "HashPtrTable.H"
#include "motionSolver.H"
#include "tetMotionSolver.H"
#include "fvMotionSolver.H"
#include "regIOobject.H"
#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "movingInterfacePatches.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class interfaceTrackingFvMesh Declaration
\*---------------------------------------------------------------------------*/

class interfaceTrackingFvMesh
:
    virtual public topoChangerFvMesh
{
    // Private Data

        // Dynamic motion solver
        autoPtr<motionSolver> motionPtr_;

        // Motion Solver type
        bool fvMotionSolver_;
        bool feMotionSolver_;

        //- Motion dictionary
        dictionary motionDict_;

        //- Moving interface patches
        HashPtrTable<movingInterfacePatches> movingInterfaces_;

        PtrList<entry> movingInterfaceEntries_;

        //- Moving / non-moving
        bool active_;


    // Private Member Functions

        //- No copy construct
        interfaceTrackingFvMesh(const interfaceTrackingFvMesh&) = delete;

        //- No copy assignment
        void operator=(const interfaceTrackingFvMesh&) = delete;

public:

    //- Runtime type information
    TypeName("interfaceTrackingFvMesh");


    // Constructors

        //- Construct from IOobject
        interfaceTrackingFvMesh(const IOobject& io);


    //- Destructor
    ~interfaceTrackingFvMesh();


    // Member Functions

        //- Return access to dictionary
        inline const dictionary& motionDict() const
        {
            return motionDict_;
        }


        //- Return reference to finite volume mesh
        fvMesh& mesh()
        {
            return *this;
        }

        //- Return reference to finite volume mesh
        const fvMesh& mesh() const
        {
            return *this;
        }

        bool& active()
        {
            return active_;
        }

        //- Return constant reference to velocity field
        const volVectorField& U() const;

        //- Return constant reference to pressure field
        const volScalarField& p() const;

        //- Return constant reference to velocity field
        const surfaceScalarField& phi() const;

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
