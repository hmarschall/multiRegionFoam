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

Class
    regionInterface

Description
    Class to handle the coupling of region patches.

SourceFiles
    regionInterface.C

\*---------------------------------------------------------------------------*/

#ifndef regionInterface_H
#define regionInterface_H

#include "Ostream.H"
#include "IOdictionary.H"
#include "dynamicFvMesh.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "foamTime.H"

#include "interfaceKey.H"

//#include "volFields.H"
//#include "primitivePatchInterpolation.H"
//#include "fvPatch.H"
#include "standAlonePatch.H"
//#include "fvPatchFields.H"
//#include "volFields.H"
//#include "fvPatchFieldMapper.H"
#include "fvCFD.H"
#include "faCFD.H"

#include "globalPolyPatch.H"
#include "regionCouplePolyPatch.H"

#include "ggiInterpolation.H"
#include "RBFInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class regionInterface Declaration
\*---------------------------------------------------------------------------*/

class regionInterface
:
    public IOdictionary,
    public interfaceKey
{
    // Private data

        //- multiRegion properties dictionary
        IOdictionary multiRegionProperties_;

        //- interface properties dictionary
        IOdictionary regionInterfaceProperties_;

        //- Reference to time
        const Time& runTime_;

        //- Reference to the fvPatches
        const fvPatch& patchA_;

        const fvPatch& patchB_;

        //- Reference to fvMesh A
        const fvMesh& meshA_;

        //- Reference to fvMesh B
        const fvMesh& meshB_;

        //- Swicth : attached/detached state of mesh interface
        Switch attachedA_;

        Switch attachedB_;

        //- Interpolator update frequency
        label interpolatorUpdateFrequency_;

        //- A face zone current points
        mutable autoPtr<vectorField> currentAZonePointsPtr_;

        //- A face zone primitive patch in current configuration
        mutable autoPtr<standAlonePatch> currentAZonePatchPtr_;

        //- GGI interpolator between two standAlone patches
        mutable autoPtr< GGIInterpolation<standAlonePatch, standAlonePatch> >
            ggiInterpolatorPtr_;

        //- Global poly patch
        //  This is a copy of a patch which is reconstructed in full on each
        //  processor. This is convenient for interpolating between two patches
        //  in parallel
        //- A-side
        mutable autoPtr<globalPolyPatch> globalPatchAPtr_;
        //- B-side
        mutable autoPtr<globalPolyPatch> globalPatchBPtr_;
        
        //- finite area mesh
        mutable autoPtr<faMesh> aMeshPtr_;

protected:

        //- surface velocity field
        mutable autoPtr<areaVectorField> UsPtr_;

        //- interfacial curvature
        mutable autoPtr<areaScalarField> KPtr_;

        //- surface flux field
        mutable autoPtr<edgeScalarField> phisPtr_;


private:

    // Private Member Functions

        //- Calculate A-to-B / B-to-A ggi interpolator
        void calcGgiInterpolator() const;

        //- Calculate current A zone primitive patch
        void calcCurrentAZonePatch() const;

        //- Calculate current A zone points
        void calcCurrentAZonePoints() const;

        //- Disallow default bitwise copy construct
        regionInterface(const regionInterface&);

        //- Disallow default bitwise assignment
        void operator=(const regionInterface&);

public:

    //- Runtime type information
    TypeName("regionInterface");


    // Declare run-time constructor selection table
#ifndef SWIG
        declareRunTimeSelectionTable
        (
            autoPtr,
            regionInterface,
            IOdictionary,
            (
                const Time& runTime,
                const fvPatch& patchA,
                const fvPatch& patchB
            ),
            (runTime, patchA, patchB)
        );
#endif

    // Constructors

        //- Construct from components
        regionInterface
        (
//            const word& type,
            const Time& runTime,
            const fvPatch& patchA,
            const fvPatch& patchB
        );


    // Selectors

        //- Select constructed from fluid and solid meshes
        static autoPtr<regionInterface> New
        (
            const Time& runTime,
            const fvPatch& patchA,
            const fvPatch& patchB
        );


    // Destructor

        virtual ~regionInterface(){}


    // Public enumerations

        //- Method for transferring information between the solid and the fluid
        //  interfaces
        enum interfaceTransferMethod
        {
            directMap,
            RBF,
            GGI
        };

        //- interfaceTransferMethod names
        static const NamedEnum<interfaceTransferMethod, 3>
            interfaceTransferMethodNames_;

private:
        //- Transfer method
        interfaceTransferMethod transferMethod_;

public:
    // Member Functions

        // Access

            //- Return reference to region properties
            const IOdictionary& regionProperties()
            {
                return multiRegionProperties_;
            }

            //- Return reference to interface properties
            const IOdictionary& interfaceProperties()
            {
                return regionInterfaceProperties_;
            }

			//- Return reference to DB
			const Time& runTime() const
			{
				return runTime_;
			}

            //- Return interface name
            virtual word name() const;

            //- Return a reference to the patch A
            inline const fvPatch& patchA() const
            {
                return patchA_;
            };

            //- Return a reference to the patch B
            inline const fvPatch& patchB() const
            {
                return patchB_;
            };

			//- Return reference to fvMesh es
			const fvMesh& meshA() const
			{
				return meshA_;
			}

            label patchAID() const
            {
                return patchA_.index();
            }

			const fvMesh& meshB() const
			{
				return meshB_;
			}

            label patchBID() const
            {
                return patchB_.index();
            }

            //- Return reference to faMesh
            faMesh& aMesh()
            {
                if (aMeshPtr_.empty())
                {
                    makeFaMesh();
                }
                
                return aMeshPtr_();
            }

            const faMesh& aMesh() const
            {
                if (aMeshPtr_.empty())
                {
                    makeFaMesh();
                }
                
                return aMeshPtr_();
            }
        
            // return surface velocity    
            Foam::areaVectorField& Us() 
            {
                if (UsPtr_.empty())
                {
                    makeUs();
                }

                return UsPtr_();
            }

            const Foam::areaVectorField& Us() const
            {
                if (UsPtr_.empty())
                {
                    makeUs();
                }

                return UsPtr_();
            }

            // return surface flux
            edgeScalarField& Phis()
            {
                if (phisPtr_.empty())
                {
                    makePhis();
                }

                return phisPtr_();
            }

            const edgeScalarField& Phis() const
            {
                if (phisPtr_.empty())
                {
                   makePhis();
                }

                return phisPtr_();
            }

            // return interfacial curvature
            Foam::areaScalarField& K()
            {
                if (KPtr_.empty())
                {
                    makeK();
                }

                return KPtr_();
            }

            const Foam::areaScalarField& K() const
            {
                if (KPtr_.empty())
                {
                    makeK();
                }

                return KPtr_();
            }

            //- Return GGI interpolator (B-to-A)
            const GGIInterpolation
            <
                standAlonePatch, standAlonePatch
            >& ggiInterpolator() const;

            //- Make the global poly patch
            void makeGlobalPatches() const;

            //- Return const-access to the global poly patch A
            const globalPolyPatch& globalPatchA() const;

            //- Return const-access to the global poly patch B
            const globalPolyPatch& globalPatchB() const;

            //- Clear global poly patch
            void clearGlobalPatches() const;

            //- Return current A zone patch
            const standAlonePatch& currentAZonePatch() const;

            //- Return current A zone points
            const vectorField& currentAZonePoints() const;

            //- Transfer/map/interpolate from one zone faces to another zone
            //  faces
            template<class Type>
            void transferFacesZoneToZone
            (
                const word& fromRegion,          // from region name
                const word& toRegion,            // to region name
                const standAlonePatch& fromZone, // from zone
                const standAlonePatch& toZone,   // to zone
                const Field<Type>& fromField,    // from field
                Field<Type>& toField             // to field
            ) const;

            //- Transfer/map/interpolate from zone A faces to zone B
            //  faces and vice versa
            template<class Type>
            void transferFacesAToB
            (
                const Field<Type>& fromField,    // from field
                Field<Type>& toField             // to field
            ) const;

            template<class Type>
            void transferFacesBToA
            (
                const Field<Type>& fromField,    // from field
                Field<Type>& toField             // to field
            ) const;

            template<class Type>
            tmp<Field<Type> > transferFacesFromA
            (
                const Field<Type>& fromField
            ) const;

        // Edit

            //- Update interpolator
            virtual void updateInterpolatorAndGlobalPatches();

            //- Attach the mesh interface using mesh modifier
            virtual void attach();

            //- Detach the mesh interface using mesh modifier
            virtual void detach();

            // interface-specific physics
            void makeFaMesh() const;

            virtual void makeUs() const = 0;
            virtual void makeK() const = 0;
            virtual void makePhis() const = 0;

        // Write

           //- Write entries for patches
//           void writeEntries(Ostream& os) const;
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "regionInterfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
