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

#include "fixedGradientFaPatchFields.H"
#include "regionInterfaceType.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionInterfaceType, 0);
    defineRunTimeSelectionTable(regionInterfaceType, IOdictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

word Foam::regionInterfaceType::assembleName
(
    const fvPatch& patchA,
    const fvPatch& patchB,
    const word& typeName
)
{
    word meshAName = patchA.boundaryMesh().mesh().name();

    word patchAName = patchA.name();
    word PatchAName = word(toupper(patchAName[0]));
    PatchAName = PatchAName + patchAName.substr(1);

    word meshBName = patchB.boundaryMesh().mesh().name();
    word MeshBName = word(toupper(meshBName[0]));
    MeshBName = MeshBName + meshBName.substr(1);

    word patchBName = patchB.name();
    word PatchBName = word(toupper(patchBName[0]));
    PatchBName = PatchBName + patchBName.substr(1);

    word InterfaceTypeName = word(toupper(typeName[0]));
    InterfaceTypeName = InterfaceTypeName + typeName.substr(1);

    return
    (
        meshAName + PatchAName
      + MeshBName + PatchBName
//      + InterfaceTypeName
    );
}

void Foam::regionInterfaceType::makeGlobalPatches() const
{
    if (globalPatchAPtr_.valid() || globalPatchBPtr_.valid())
    {
        FatalErrorIn(type() + "::makeGlobalPatches() const")
            << "Pointer already set!" << abort(FatalError);
    }

    Info<< "Creating global patches : "
    << patchA().name() << " and "
    << patchB().name() << " for regionInterfaceType "
    << interfaceName()
    << endl;

    globalPatchAPtr_.set(new globalPolyPatch(patchA().name(), meshA()));
    globalPatchBPtr_.set(new globalPolyPatch(patchB().name(), meshB()));

    globalPatchAPtr_().globalPatch();
    globalPatchBPtr_().globalPatch();
}

void Foam::regionInterfaceType::clearGlobalPatches() const
{
    globalPatchAPtr_.clear();
    globalPatchBPtr_.clear();
}

void Foam::regionInterfaceType::makeInterfaceToInterface() const
{
    if (interfaceToInterfacePtr_.valid())
    {
        FatalErrorIn
        (
            "void Foam::interfaceToInterfaceMapping::"
            "makeInterfaceToInterface() const"
        )   << "Mapping object already set!" << abort(FatalError);
    }

    // Lookup the type
    const word type = regionInterfaceProperties_.lookupOrDefault<word>
    (
        "interfaceTransferMethod", "GGI"
    );

    interfaceToInterfacePtr_ =
    (
        interfaceToInterfaceMapping::New
        (
            type,
            regionInterfaceProperties_.subDict(type + "Coeffs"),
            meshA().boundaryMesh()[patchAID()],
            meshB().boundaryMesh()[patchBID()],
            globalPatchA(),
            globalPatchB()
        )
    );
}

void Foam::regionInterfaceType::resetFaMesh() const
{
    word aMeshName = patchA().name() + "FaMesh";

    // look up faMesh from object registry
    if (meshA().foundObject<faMesh>(aMeshName))
    {
        aMeshPtr_.reset
        (
            const_cast<faMesh*>
            (
                &meshA().lookupObject<faMesh>(aMeshName)
            )
        );
    }
}

void Foam::regionInterfaceType::makeFaMesh() const
{
    if (!aMeshPtr_.empty())
    {
        FatalErrorIn("regionInterfaceType::makeFaMesh()")
            << "finite area mesh already exists"
            << abort(FatalError);
    }

    word aMeshName = patchA().name() + "FaMesh";

    // look up faMesh from object registry
    if (meshA().foundObject<faMesh>(aMeshName))
    {
        aMeshPtr_.reset
        (
            const_cast<faMesh*>
            (
                &meshA().lookupObject<faMesh>(aMeshName)
            )
        );
    }
    // or create new mesh
    else
    {
        aMeshPtr_.set(new faMesh(meshA()));

        // set unique name
        aMeshPtr_->rename(patchA().name() + "FaMesh");
    }
}

void Foam::regionInterfaceType::correctCurvature
(
    areaScalarField& K
)
{
    scalarField& KI = K.internalField();

    forAll(curvatureCorrectedSurfacePatches_, patchI)
    {
        label patchID =
            aMesh().boundary().findPatchID
            (
                curvatureCorrectedSurfacePatches_[patchI]
            );

        if(patchID == -1)
        {
            FatalErrorIn("regionInterfaceType::correctCurvature(...)")
                << "Wrong faPatch name in the curvatureCorrectedSurfacePatches"
                    << " list defined in regionInterfaceProperties"
                    << abort(FatalError);
        }

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

void regionInterfaceType::clearOut() const
{
    interfaceToInterfacePtr_.clear();
//    aMeshPtr_.clear();

    clearGlobalPatches();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionInterfaceType::regionInterfaceType
(
    const word& type,
    const dictionary& dict,
    const Time& runTime,
    const fvPatch& patchA,
    const fvPatch& patchB
)
:
    IOdictionary
    (
        IOobject
        (
            assembleName(patchA, patchB, type),
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    MeshObject<fvMesh, regionInterfaceType>(patchA.boundaryMesh().mesh()),
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
    globalPatchAPtr_(),
    globalPatchBPtr_(),
    interfaceToInterfacePtr_(),
    meshA_(patchA_.boundaryMesh().mesh()),
    meshB_(patchB_.boundaryMesh().mesh()),
    attachedA_(false),
    attachedB_(false),
    changing_(false),
    moving_(false),
    interpolatorUpdateFrequency_
    (
        regionInterfaceProperties_
        .lookupOrDefault<int>("interpolatorUpdateFrequency", 1)
    ),
    aMeshPtr_(), //new faMesh(meshA_)
//    areaMesh_(faMesh(meshA_)),
    curvatureCorrectedSurfacePatches_
    (
        regionInterfaceProperties_.lookup("curvatureCorrectedSurfacePatches")
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

    // Force creation of interface-to-interface object
    // as they may need to read fields on restart
    interfaceToInterface();

    if (debug)
    {
        //Output region interface information
        Pout<< "regionInterfaceType Info: " << interfaceName() << nl
            << "local patchA: " << nl
            << " name: " << patchA_.name()
            << " size: " << patchA_.size()
            << " nPoints: " << patchA_.patch().nPoints()
            << " nEdges: " << patchA_.patch().nEdges()
            << nl
            << "local patchB: " << nl
            << " name: " << patchB_.name()
            << " size: " << patchB_.size()
            << " nPoints: " << patchB_.patch().nPoints()
            << " nEdges: " << patchB_.patch().nEdges()
            << nl
            << "global patchA: " << nl
            << " name: " << globalPatchAPtr_->patchName()
            << " size: " << globalPatchAPtr_->globalPatch().size()
            << " nPoints: " << globalPatchAPtr_->globalPatch().nPoints()
            << " nEdges: " << globalPatchAPtr_->globalPatch().nEdges()
            << nl
            << "global patchB: " << nl
            << " name: " << globalPatchBPtr_->patchName()
            << " size: " << globalPatchBPtr_->globalPatch().size()
            << " nPoints: " << globalPatchBPtr_->globalPatch().nPoints()
            << " nEdges: " << globalPatchBPtr_->globalPatch().nEdges()
            << nl
            << endl;
    }
}

// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

regionInterfaceType::~regionInterfaceType()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::regionInterfaceType::interfaceName() const
{
    return IOdictionary::name();
}

const Foam::globalPolyPatch& Foam::regionInterfaceType::globalPatchA() const
{
    if (globalPatchAPtr_.empty())
    {
        FatalErrorIn(type() + "::makeGlobalPatch() const")
            << "makeGlobalPatches(...) must be called before globalPatchA "
            << "can be called!" << abort(FatalError);
    }

    return globalPatchAPtr_();
}

const Foam::globalPolyPatch& Foam::regionInterfaceType::globalPatchB() const
{
    if (globalPatchBPtr_.empty())
    {
        FatalErrorIn(type() + "::makeGlobalPatch() const")
            << "makeGlobalPatches(...) must be called before globalPatchA "
            << "can be called!" << abort(FatalError);
    }

    return globalPatchBPtr_();
}

void Foam::regionInterfaceType::updateInterpolatorAndGlobalPatches()
{
    Info << "Updating interpolator and global patches for regionInterfaceType" << endl;

    if (interfaceToInterfacePtr_.empty())
    {
        interfaceToInterface();
    }
    else if (interpolatorUpdateFrequency_ != 0)
    {
        if
        (
            ((runTime().timeIndex() - 1) % interpolatorUpdateFrequency_) == 0
         || (!moving() && changing())
        )
        {
            // Clear current interpolators
            interfaceToInterfacePtr_.clear();

            // Clear and re-create global patches
            clearGlobalPatches();
            makeGlobalPatches();

            // Re-create interpolators
            interfaceToInterface();
        }
    }
}

const Foam::interfaceToInterfaceMapping&
Foam::regionInterfaceType::interfaceToInterface() const
{
    if (interfaceToInterfacePtr_.empty())
    {
        makeInterfaceToInterface();
    }

    return interfaceToInterfacePtr_();
}

void Foam::regionInterfaceType::attach()
{
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

void Foam::regionInterfaceType::detach()
{
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

void Foam::regionInterfaceType::update()
{
    // critical: only if mesh topology has changed
    if (!moving() && changing())
    {
        clearOut();
        resetFaMesh();
        makeGlobalPatches();
    }

    updateK();
    this->correct();
}

void Foam::regionInterfaceType::updateK()
{
    areaScalarField& curv =
        const_cast<areaScalarField&>
        (
           aMesh().faceCurvatures()
        );

    correctCurvature(curv);

    curv.correctBoundaryConditions();
}

// Update for mesh motion
bool Foam::regionInterfaceType::movePoints() const
{
    return true;
}


// Update on topology change
bool Foam::regionInterfaceType::updateMesh(const mapPolyMesh& mpm) const
{
    if (debug)
    {
        Info << "Clearing out regionInterfaceType after topology change" << endl;
    }

    // Wipe out demand-driven data
    clearOut();
//    resetFaMesh();
//    makeGlobalPatches();
//    interfaceToInterface();

    return true;
}

// ************************************************************************* //
