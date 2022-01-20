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

\*---------------------------------------------------------------------------*/

#include "regionInterface.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionInterface, 0);
    defineRunTimeSelectionTable(regionInterface, IOdictionary);

    template<>
    const char*
    NamedEnum<regionInterface::interfaceTransferMethod, 3>::names[] =
    {
        "directMap",
        "RBF",
        "GGI"
    };

    const NamedEnum<regionInterface::interfaceTransferMethod, 3>
        regionInterface::interfaceTransferMethodNames_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::regionInterface::calcCurrentAZonePatch() const
{
    // Find global face zones
    if (!currentAZonePatchPtr_.empty())
    {
        FatalErrorIn
        (
            "void regionInterface::calcCurrentAZonePatch() const"
        )   << "Current A zone patch already exists"
            << abort(FatalError);
    }

    currentAZonePatchPtr_.set
    (
        new standAlonePatch
        (
            globalPatchA().globalPatch().localFaces(),
            currentAZonePoints()
        )
    );
}

void Foam::regionInterface::calcCurrentAZonePoints() const
{
    // Find global face zones
    if (!currentAZonePointsPtr_.empty())
    {
        FatalErrorIn
        (
            "void regionInterface::"
            "calcCurrentSolidAPoints() const"
        )   << "Current A zone points already exist"
            << abort(FatalError);
    }

    // Calculate global patch deformed points
    tmp<Foam::vectorField> currentFaceZonePoints =
        globalPatchA().globalPatch().localPoints();

    if
    (
        meshA().objectRegistry::foundObject<vectorIOField>("totalDisplacement")
    )
    {
        // Patch point displacement
        const vectorIOField& pointDisplacement =
            meshA().lookupObject<vectorIOField>("totalDisplacement");

        currentFaceZonePoints() += 
            globalPatchA().patchPointToGlobal(pointDisplacement);
    }

//    tmp<Foam::vectorField> currentFaceZonePoints =
//        globalPatchA().globalPatch().localPoints();
//      + globalPatchA().patchPointToGlobal(pointDisplacement);

    // Return current A zone points
    currentAZonePointsPtr_.set
    (
        new vectorField(currentFaceZonePoints())
    );
}

void Foam::regionInterface::calcGgiInterpolator() const
{
    // Create ggi interpolation
    if (ggiInterpolatorPtr_.valid())
    {
        FatalErrorIn
        (
            "void regionInterface::"
            "calcGgiInterpolator() const"
        )   << "Ggi interpolator already exists"
            << abort(FatalError);
    }

    // Remove current A zone and points so that it will be re-created in the
    // deformed position
    currentAZonePatchPtr_.clear();
    currentAZonePointsPtr_.clear();

    Info<< "Create GGI zone-to-zone interpolator" << endl;

    ggiInterpolatorPtr_.set
    (
        new GGIInterpolation<standAlonePatch, standAlonePatch>
        (
            globalPatchB().globalPatch(),
            currentAZonePatch(),
            tensorField(0),
            tensorField(0),
            vectorField(0), // Slave-to-master separation. Bug fix
            true,           // Patch data is complete on all processors
            SMALL,          // Non-overlapping face tolerances
            SMALL,
            true,           // Rescale weighting factors
            ggiInterpolation::BB_OCTREE
        )
    );

    Info<< "Checking A-to-B point interpolator" << endl;
    {
        const vectorField AZonePointsAtB =
            ggiInterpolatorPtr_().slaveToMasterPointInterpolate
            (
                currentAZonePoints()
            );

        const vectorField BZonePoints =
            globalPatchB().globalPatch().localPoints();

        const scalar maxDist = gMax
        (
            mag
            (
                BZonePoints
              - AZonePointsAtB
            )
        );

        Info<< "A-to-B point interpolation error: " << maxDist
            << endl;
    }

    Info<< "Checking B-to-A face interpolator" << endl;
    {
        const vectorField BPatchFaceCentres =
            vectorField
            (
                meshB().boundaryMesh()[patchBID()].faceCentres()
            );

        const vectorField BZoneFaceCentres =
            globalPatchB().patchFaceToGlobal(BPatchFaceCentres);

        const vectorField AZoneFaceCentres =
            ggiInterpolatorPtr_().masterToSlave
            (
                BZoneFaceCentres
            );

        const vectorField APatchFaceCentres =
            globalPatchA().globalFaceToPatch(AZoneFaceCentres);

        scalar maxDist = gMax
        (
            mag
            (
                APatchFaceCentres
              - meshA().boundaryMesh()[patchA().index()].faceCentres()
            )
        );

        Info<< "B-to-A face interpolation error: " << maxDist
            << endl;
    }

    Info<< "Number of uncovered master faces: "
        << ggiInterpolatorPtr_().uncoveredMasterFaces().size() << endl;

    Info<< "Number of uncovered slave faces: "
        << ggiInterpolatorPtr_().uncoveredSlaveFaces().size() << endl;

    ggiInterpolatorPtr_().slavePointDistanceToIntersection();
    ggiInterpolatorPtr_().masterPointDistanceToIntersection();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionInterface::regionInterface
(
    const Time& runTime,
    const fvPatch& patchA,
    const fvPatch& patchB
)
:
    IOdictionary
    (
        IOobject
        (
            patchA.boundaryMesh().mesh().name() + patchA.name() 
            + patchB.boundaryMesh().mesh().name() + patchB.name(),
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    interfaceKey(patchA.name(), patchB.name()),
    multiRegionProperties_
    (
        IOobject
        (
            "multiRegionProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    regionInterfaceProperties_
    (
        IOobject
        (
            "regionInterfaceProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    gravitationalProperties_
    (
        IOobject
        (
            "g",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),    
    runTime_(runTime),
    patchA_(patchA),
    patchB_(patchB),
    meshA_(patchA_.boundaryMesh().mesh()),
    meshB_(patchB_.boundaryMesh().mesh()),
    attachedA_(false),
    attachedB_(false),
    interpolatorUpdateFrequency_
    (
        multiRegionProperties_.lookupOrDefault<int>("interpolatorUpdateFrequency", 0)
    ),
    currentAZonePointsPtr_(),
    currentAZonePatchPtr_(),
    ggiInterpolatorPtr_(),
    globalPatchAPtr_(),
    globalPatchBPtr_(),
    aMeshPtr_(),
    UsPtr_(),
    KPtr_(),
    phisPtr_(),
    transferMethod_
    (
        interfaceTransferMethodNames_
        [
            multiRegionProperties_.lookupOrDefault<word>
            (
                "interfaceTransferMethod", "GGI"
            )
        ]
    )
{
    // Create global patches
    makeGlobalPatches();

    // Get initial state of coupled patches on A/B side
    const polyPatchList& patchesA = meshA().boundaryMesh();
    const polyPatchList& patchesB = meshB().boundaryMesh();

    if (isType<regionCouplePolyPatch>(patchesA[patchAID()]))
    {
        const regionCouplePolyPatch& rcp =
            refCast<const regionCouplePolyPatch>(patchesA[patchAID()]);

            // Check if coupled
            if (rcp.coupled())
            {
                attachedA_ = true;
            }
    }

    if (isType<regionCouplePolyPatch>(patchesB[patchBID()]))
    {
        const regionCouplePolyPatch& rcp =
            refCast<const regionCouplePolyPatch>(patchesB[patchBID()]);

            // Check if coupled
            if (rcp.coupled())
            {
                attachedB_ = true;
            }
    }

    Info << "This is the regionInterface : " << name() << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::regionInterface::name() const
{
    word meshAName = meshA_.name();

    word meshBName = meshB_.name();
    meshBName[0] = toupper(meshBName[0]);

    word name1(Pair<word>::first());
    name1[0] = toupper(name1[0]);

    word name2(Pair<word>::second());
    name2[0] = toupper(name2[0]);

    return meshAName + name1 + meshBName + name2;
}

void Foam::regionInterface::makeGlobalPatches() const
{
    if (globalPatchAPtr_.valid() || globalPatchBPtr_.valid())
    {
        FatalErrorIn(type() + "::makeGlobalPatches() const")
            << "Pointer already set!" << abort(FatalError);
    }

    globalPatchAPtr_.set(new globalPolyPatch(patchA().name(), meshA()));
    globalPatchBPtr_.set(new globalPolyPatch(patchB().name(), meshB()));
}

const Foam::globalPolyPatch& Foam::regionInterface::globalPatchA() const
{
    if (globalPatchAPtr_.empty())
    {
        FatalErrorIn(type() + "::makeGlobalPatch() const")
            << "makeGlobalPatches(...) must be called before globalPatchA "
            << "can be called!" << abort(FatalError);
    }

    return globalPatchAPtr_();
}

const Foam::globalPolyPatch& Foam::regionInterface::globalPatchB() const
{
    if (globalPatchBPtr_.empty())
    {
        FatalErrorIn(type() + "::makeGlobalPatch() const")
            << "makeGlobalPatches(...) must be called before globalPatchA "
            << "can be called!" << abort(FatalError);
    }

    return globalPatchBPtr_();
}

void Foam::regionInterface::clearGlobalPatches() const
{
    globalPatchAPtr_.clear();
    globalPatchBPtr_.clear();
}

const Foam::standAlonePatch&
Foam::regionInterface::currentAZonePatch() const
{
    if (currentAZonePatchPtr_.empty())
    {
        calcCurrentAZonePatch();
    }

    return currentAZonePatchPtr_();
}

const Foam::vectorField&
Foam::regionInterface::currentAZonePoints() const
{
    if (currentAZonePointsPtr_.empty())
    {
        calcCurrentAZonePoints();
    }

    return currentAZonePointsPtr_();
}

void Foam::regionInterface::updateInterpolatorAndGlobalPatches()
{
    Info << "Updating interpolator and global patches" << endl;

    if (ggiInterpolatorPtr_.empty())
    {
        ggiInterpolator();
    }
    else if (interpolatorUpdateFrequency_ != 0)
    {
        if (((runTime().timeIndex() - 1) % interpolatorUpdateFrequency_) == 0)
        {
//            deleteDemandDrivenData(ggiInterpolatorPtr_);
            clearGlobalPatches();
            makeGlobalPatches();
            ggiInterpolator();
        }
    }
}

void Foam::regionInterface::attach()
{
//    Info << "attach() : " << "meshA = " << meshA().name() << endl;
//    Info << "attach() : " << "meshB = " << meshB().name() << endl;

//    if (attachedA_ && attachedB_)
//    {
//        FatalErrorIn
//        (
//            "void Foam::regionInterface::attach()"
//        )   << "Attempt to attach. Patches already in attached mode."
//            << abort(FatalError);
//    }

    if (!attachedA_)
    {
        fvMesh& mesh = const_cast<fvMesh&>(meshA());

        const polyPatchList& patches = mesh.boundaryMesh();

        if (isType<regionCouplePolyPatch>(patches[patchAID()]))
        {
            const regionCouplePolyPatch& rcp =
                refCast<const regionCouplePolyPatch>(patches[patchAID()]);

            // Attach it here
            attachedA_ = true;
            if (rcp.master()) rcp.attach();
        }

        // Force recalculation of weights
        mesh.surfaceInterpolation::movePoints();
    }

    if (!attachedB_)
    {
        fvMesh& mesh = const_cast<fvMesh&>(meshB());

        const polyPatchList& patches = mesh.boundaryMesh();

        if (isType<regionCouplePolyPatch>(patches[patchBID()]))
        {
            const regionCouplePolyPatch& rcp =
                refCast<const regionCouplePolyPatch>(patches[patchBID()]);

            // Attach it here
            attachedB_ = true;
            if (rcp.master()) rcp.attach();
        }

        // Force recalculation of weights
        mesh.surfaceInterpolation::movePoints();
    }
}

void Foam::regionInterface::detach()
{
//    Info << "detach() : " << "meshA = " << meshA().name() << endl;
//    Info << " patch name A = " << patchAName() << endl;

//    Info << "detach() : " << "meshB = " << meshB().name() << endl;
//    Info << " patch name B = " << patchBName() << endl;


//    if (!attachedA_ && !attachedB_)
//    {
//        FatalErrorIn
//        (
//            "void Foam::regionInterface::detach()"
//        )   << "Attempt to detach. Patches already in detached mode."
//            << abort(FatalError);
//    }

    if (attachedA_)
    {
        fvMesh& mesh = const_cast<fvMesh&>(meshA());

        const polyPatchList& patches = mesh.boundaryMesh();

        if (isType<regionCouplePolyPatch>(patches[patchAID()]))
        {
            const regionCouplePolyPatch& rcp =
                refCast<const regionCouplePolyPatch>(patches[patchAID()]);

            // Detach it here
            attachedA_ = false;
            if (rcp.master()) rcp.detach();
        }

        // Force recalculation of weights
        mesh.surfaceInterpolation::movePoints();
    }

    if (attachedB_)
    {
        fvMesh& mesh = const_cast<fvMesh&>(meshB());

        const polyPatchList& patches = mesh.boundaryMesh();

        if (isType<regionCouplePolyPatch>(patches[patchBID()]))
        {
            const regionCouplePolyPatch& rcp =
                refCast<const regionCouplePolyPatch>(patches[patchBID()]);

            // Detach it here
            attachedB_ = false;
            if (rcp.master()) rcp.detach();
        }

        // Force recalculation of weights
        mesh.surfaceInterpolation::movePoints();
    }
}

const Foam::GGIInterpolation<standAlonePatch, standAlonePatch>&
Foam::regionInterface::ggiInterpolator() const
{
    if (ggiInterpolatorPtr_.empty())
    {
        calcGgiInterpolator();
    }

    return ggiInterpolatorPtr_();
}

void Foam::regionInterface::makeFaMesh() const
{
    if (!aMeshPtr_.empty())
    {
        FatalErrorIn("regionInterface::makeFaMesh()")
            << "finite area mesh already exists"
            << abort(FatalError);
    }

    aMeshPtr_.set(new faMesh(meshA())); 
}  

void Foam::regionInterface::makeK() const
{
    if (!KPtr_.empty())
    {
        FatalErrorIn("regionInterface::makeK()")
            << "surface curvature field already exists"
            << abort(FatalError);
    }

    KPtr_.set
    (
        new areaScalarField
        (
            IOobject
            (
                "K",
                runTime().constant(), 
                runTime(), 
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh(),
            dimensioned<scalar>("K", dimless/dimLength, pTraits<scalar>::zero),
            zeroGradientFaPatchVectorField::typeName
        )
    );
}


void Foam::regionInterface::makePhis() const
{

    const surfaceScalarField& phi = 
        meshA().lookupObject<surfaceScalarField>("phi");

    if (!phisPtr_.empty())
    {
        FatalErrorIn("regionInterface::makePhis()")
            << "surface fluid flux already exists"
            << abort(FatalError);
    }

    phisPtr_.set
    (
        new edgeScalarField
        (
            IOobject
            (
                phi.name() + "s",
                runTime().constant(), 
                runTime(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            linearEdgeInterpolate(Us()) & aMesh().Le()
        )
    );
} 


void Foam::regionInterface::makeUs() const
{
    // error if U is initialized in the constructor
    const volVectorField& U = meshA().lookupObject<volVectorField>("U");
         
    if (!UsPtr_.empty())
    {
        FatalErrorIn("regionInterface::makeUs()")
            << "surface velocity field already exists"
            << abort(FatalError);
    }

    wordList patchFieldTypes
    (
        aMesh().boundary().size(),
        zeroGradientFaPatchVectorField::typeName
    );

    forAll(aMesh().boundary(), patchI) 
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == wedgeFaPatch::typeName
        )
        {
            patchFieldTypes[patchI] = 
                wedgeFaPatchVectorField::typeName;
        }
        else
        {
            label ngbPolyPatchID = 
                aMesh().boundary()[patchI].ngbPolyPatchIndex();

            if (ngbPolyPatchID != -1)
            {
                if
                (
                    meshA().boundary()[ngbPolyPatchID].type() 
                 == wallFvPatch::typeName
                )
                {
                    patchFieldTypes[patchI] =
                        slipFaPatchVectorField::typeName;
                }
            }
        }
    }
    
    UsPtr_.set
    (
        new areaVectorField
        (
            IOobject
            (
                U.name() + "s",
                runTime().constant(), 
                runTime(), 
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh(),
            dimensioned<vector>("Us", dimVelocity, vector::zero),
            patchFieldTypes
        )
    );
}  

const vectorField& Foam::regionInterface::Up()
{
    const volVectorField& U = meshA().lookupObject<volVectorField>("U");
        
    const fvBoundaryMesh& fvbm = meshA().boundary(); 

    const fvPatch& p = fvbm[patchAID()];

    return p.lookupPatchField<volVectorField, vector>(U.name());
} 


void Foam::regionInterface::correctUsBoundaryConditions()
{  
    const volVectorField& U = meshA().lookupObject<volVectorField>("U");
         
    forAll(Us().boundaryField(), patchI)
    {
        if
        (
            UsPtr_().boundaryField()[patchI].type()
         == calculatedFaPatchVectorField::typeName
        )
        {
            vectorField& pUs = Us().boundaryField()[patchI];

            pUs = Us().boundaryField()[patchI].patchInternalField();

            label ngbPolyPatchID =
                aMesh().boundary()[patchI].ngbPolyPatchIndex();

            if (ngbPolyPatchID != -1)
            {
                if
                (
                    (
                        U.boundaryField()[ngbPolyPatchID].type()
                     == slipFvPatchVectorField::typeName
                    )
                 ||
                    (
                        U.boundaryField()[ngbPolyPatchID].type()
                     == symmetryFvPatchVectorField::typeName
                    )
                )
                {
                    vectorField N
                    (
                        aMesh().boundary()[patchI].ngbPolyPatchFaceNormals()
                    );

                    pUs -= N*(N&pUs);
                }
            }
        }
    }

    Us().correctBoundaryConditions();
}


void Foam::regionInterface::updateUs()
{
    Us().internalField() = Up();

    correctUsBoundaryConditions();
}


void Foam::regionInterface::updateSurfaceFlux()
{
    Phis() = fac::interpolate(Us()) & aMesh().Le();
}


void Foam::regionInterface::calcCurvatureAxis
(
    areaScalarField& K
)
{
    scalarField& KI = K.internalField();

    //  TODO: dictionary lookup (also hard coded in surfaceTracking!)
    //        + use also for fixedSurfacesPatches
    label patchID = aMesh().boundary().findPatchID("centerline");

    if (patchID != -1)
    {
        const labelList& eFaces =
            aMesh().boundary()[patchID].edgeFaces();

        const labelListList& fFaces = aMesh().patch().faceFaces();

        forAll(eFaces, edgeI)
        {
            const label& curFace = eFaces[edgeI];
            const labelList& curFaceFaces = fFaces[curFace];

            scalar avrK = 0.0;
            label counter = 0;

            forAll(curFaceFaces, faceI)
            {
                label index = findIndex(eFaces, curFaceFaces[faceI]);

                if (index == -1)
                {
                    avrK += K[curFaceFaces[faceI]];
                    counter++;
                }
            }
            avrK /= counter;

            KI[curFace] = avrK;
        }
//        label counter = 0;
//        do
//        {
//            counter++;

//            K.correctBoundaryConditions();
//            areaVectorField gradK = fac::grad(K);
//            vectorField& gradKI = gradK.internalField();

//            const labelList& eFaces =
//                aMesh().boundary()[patchID].edgeFaces();

//            const labelListList& fFaces = aMesh().patch().faceFaces();

//            const vectorField& fCentres = aMesh().areaCentres();

//            forAll(eFaces, edgeI)
//            {
//                const label& curFace = eFaces[edgeI];
//                const labelList& curFaceFaces = fFaces[curFace];

//                scalar avrK = 0.0;
//                label counter = 0;

//                forAll(curFaceFaces, faceI)
//                {
//                    label index = findIndex(eFaces, curFaceFaces[faceI]);

//                    if (index == -1)
//                    {
//                        vector dr = 
//                            fCentres[curFace] 
//                          - fCentres[curFaceFaces[faceI]];

//                        avrK += KI[curFaceFaces[faceI]]
//                             + (dr&gradKI[curFaceFaces[faceI]]);
//                        counter++;
//                    }
//                }

//                avrK /= counter;

//                KI[curFace] = avrK;
//            }
//        }
//        while(counter<10);
    }
}


void Foam::regionInterface::updateK()
{
    K().internalField() = 
        const_cast<areaScalarField&>
        (
           aMesh().faceCurvatures()
        );

    calcCurvatureAxis(K());

    K().correctBoundaryConditions();
}



//void Foam::regionInterface::writeEntries(Ostream& os) const
//{
//    os.writeKeyword("regionBName");
//    os << regionBName_ << token::END_STATEMENT << nl;
//    os.writeKeyword("patchBName");
//    os << patchBName_ << token::END_STATEMENT << nl;
//    os.writeKeyword("BFieldName");
//    os << BFieldName_ << token::END_STATEMENT << nl;
//}


// ************************************************************************* //
