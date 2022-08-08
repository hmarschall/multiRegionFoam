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
    Foam::translationalMRFZone

Description
    Translational MRF zone definition based on parameters
    obtained from a control dictionary constructed from the given stream.

    The translation of the MRF region is defined by an center of
    reference frame, its velocity and acceleration.

SourceFiles
    translationalMRFZone.C

\*---------------------------------------------------------------------------*/

#ifndef translationalMRFZone_H
#define translationalMRFZone_H

#include "dictionary.H"
#include "IOReferencer.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedVector.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "fvMatrices.H"
#include "tmp.H"

#include "fixedValueFvPatchFields.H"
#include "noSlipWallFvPatchVectorField.H"
#include "mixedFvPatchFields.H"
#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class translationalMRFZone Declaration
\*---------------------------------------------------------------------------*/

class translationalMRFZone
{
    // Private data

        //- Name of the MRF zone
        const word name_;

        //- Mesh reference
        const fvMesh& mesh_;

        //- Velocity in MRF
        volVectorField Umrf_;

        //- Reference to dictionary
        const dictionary dict_;

        //- PID-controller related
        dimensionedScalar KP_;
        dimensionedScalar KI_;
        dimensionedScalar KD_;

        dimensionedVector PIDOut_;

        //- translational MRF related
        word wallPatchName_;
        word inletPatchName_;
        word spacePatchName_;
        label wallPatchID_;
        label inletPatchID_;
        label spacePatchID_;
        dimensionedVector centreCurrent_;
        dimensionedVector centreTarget_;
        dimensionedVector errorValue_;
        dimensionedVector integralComponent_;

        //- current particle velocity (relative to IRF)
        dimensionedVector UP_;

        //- current particle velocity (relative to MRF)
        dimensionedVector UPrel_;

        //- current adjustment velocity to keep particle in fixed position
        dimensionedVector dUF_;

        //- mask for adjustment velocity based on solution directions
        vector dUFMask_;

        //- mask correct centreTarget and centreCurrent 
        //  for wedges/symmetryPlanes
        vector centreMask_;

        //- current particle acceleration / relative to IRF
        dimensionedVector aP_;
        dimensionedVector aPrel_;

        //- current MRF acceleration (compared to IRF)
        dimensionedVector aF_;

        //- current velocity of the MRF (relative to IRF)
        dimensionedVector UF_;

        //- current position of the MRF (relative to IRF)
        dimensionedVector xF_;

        // moving wall velocity due to MRF
        dimensionedVector Uwall_;

        // methods to compute frame centre from
        word regionName_;
        word fieldName_;
        Switch normalise01_;
        word modelName_;

        Switch projectUWall_;
        Switch projectUInlet_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        translationalMRFZone(const translationalMRFZone&);

        //- Disallow default bitwise assignment
        void operator=(const translationalMRFZone&);

        //- Read MRF centre control
        void readMRFCentreControl();

        //- Set MRF centre control
        void setMRFCentreControl();

        //- Set wall/inlet/space patch IDs
        void setPatchIDs();

        //- Calculate PID error
        void PIDCorr();

public:

    // Declare name of the class and its debug switch
    ClassName("translationalMRFZone");

    //- Centre control options
    enum centreControls
    {
        Mesh,
        Field,
        Model
    };

    static const NamedEnum<centreControls, 3> centreControlNames_;
    centreControls centreControl_;

    // Constructors

        //- Construct from fvMesh and Istream
        translationalMRFZone(const fvMesh& mesh, Istream& is);

        //- Return clone
        autoPtr<translationalMRFZone> clone() const
        {
            notImplemented("autoPtr<translationalMRFZone> clone() const");
            return autoPtr<translationalMRFZone>(nullptr);
        }

        //- Return a pointer to a new translationalMRFZone created on freestore
        //  from Istream
        class iNew
        {
            const fvMesh& mesh_;

        public:

            iNew(const fvMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<translationalMRFZone> operator()(Istream& is) const
            {
                return autoPtr<translationalMRFZone>
                (
                    new translationalMRFZone(mesh_, is)
                );
            }
        };

    //- Destructor
    ~translationalMRFZone();


    // Member Functions

        // Access Functions

            //- Return the translationalMRFZone name
            const word& name() const
            {
                return name_;
            }

            //- Return the finite volume mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return the translationalMRFZone dictionary
            const dictionary& dict() const
            {
                return dict_;
            }

        // Used in a top level solver

            // Correct MRF properties
            void correctMRF();

            //- Add the frame acceleration contribution to the momentum equation
            void addFrameAcceleration
            (
                fvVectorMatrix& UEqn,
                const volScalarField& rho
            );

            void addFrameAcceleration
            (
                fvVectorMatrix& UEqn
            );

            //- Correct the boundary velocity for the roation of the MRF region
            void correctBoundaryVelocity
            (
                volVectorField& U,
                surfaceScalarField& phi
            );

            //- Write restart file
            void writeRestart();

            //- Read restart file
            void readRestart();

    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const translationalMRFZone&)
        {
            notImplemented
            (
                "Ostream& operator<<(Ostream& os, const translationalMRFZone&)"
            );
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
