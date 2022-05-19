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
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::regionInterface::makeGlobalPatches() const
{
    if (globalPatchAPtr_.valid() || globalPatchBPtr_.valid())
    {
        FatalErrorIn(type() + "::makeGlobalPatches() const")
            << "Pointer already set!" << abort(FatalError);
    }

    Info<< "Creating global patches : "
    << patchA().name() << " and "
    << patchB().name() << " for regionInterface "
    << name()
    << endl;

    globalPatchAPtr_.set(new globalPolyPatch(patchA().name(), meshA()));
    globalPatchBPtr_.set(new globalPolyPatch(patchB().name(), meshB()));
}

void Foam::regionInterface::clearGlobalPatches()
{
    globalPatchAPtr_.clear();
    globalPatchBPtr_.clear();
}

void Foam::regionInterface::makeInterfaceToInterface() const
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
                    WarningIn("regionInterface::makeUs() const")
                        << "Patch neighbouring to interface is wall" << nl
                        << "Not appropriate for inlets/outlets" << nl
                        << endl;

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
                runTime().timeName(), 
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
                runTime().timeName(), 
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
                runTime().timeName(), 
                runTime(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            linearEdgeInterpolate(Us()) & aMesh().Le()
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

void Foam::regionInterface::correctCurvature
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
            FatalErrorIn("regionInterface::correctCurvature(...)")
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
    globalPatchAPtr_(),
    globalPatchBPtr_(),
    interfaceToInterfacePtr_(),
    meshA_(patchA_.boundaryMesh().mesh()),
    meshB_(patchB_.boundaryMesh().mesh()),
    attachedA_(false),
    attachedB_(false),
    changing_(false),
    interpolatorUpdateFrequency_
    (
        regionInterfaceProperties_
        .lookupOrDefault<int>("interpolatorUpdateFrequency", 1)
    ),
    aMeshPtr_(),
    curvatureCorrectedSurfacePatches_
    (
        regionInterfaceProperties_.lookup("curvatureCorrectedSurfacePatches")
    ),
    UsPtr_(),
    KPtr_(),
    phisPtr_()
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

    // Force creation of interface-to-interface object 
    // as they may need to read fields on restart
    interfaceToInterface();

    if (debug)
    {
        //Output region interface information
        Pout << "regionInterface Info: " << name() << nl 
            << "    local patchA: " << nl 
            << "        name: " << patchA_.name() << " size: " << patchA_.size()  << " nPoints: " << patchA_.patch().nPoints() << " nEdges: " << patchA_.patch().nEdges()<< nl
            << "    local patchB: " << nl
            << "        name: " << patchB_.name() << " size: " << patchB_.size() << " nPoints: " << patchB_.patch().nPoints() << " nEdges: " << patchB_.patch().nEdges()<< nl  
            << "    global patchA: " << nl 
            << "        name: " << globalPatchAPtr_->patchName() << " size: " << globalPatchAPtr_->globalPatch().size()
                    << " nPoints: " << globalPatchAPtr_->globalPatch().nPoints() << " nEdges: " << globalPatchAPtr_->globalPatch().nEdges()<< nl
            << "    global patchB: " << nl 
            << "        name: " << globalPatchBPtr_->patchName() << " size: " << globalPatchBPtr_->globalPatch().size()
                    << " nPoints: " << globalPatchBPtr_->globalPatch().nPoints() << " nEdges: " << globalPatchBPtr_->globalPatch().nEdges()<< nl
            << endl;

    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::regionInterface::name() const
{
    word meshAName = meshA_.name();

    word meshBName = meshB_.name();
//    meshBName[0] = toupper(meshBName[0]);

    word name1(Pair<word>::first());
//    name1[0] = toupper(name1[0]);

    word name2(Pair<word>::second());
//    name2[0] = toupper(name2[0]);

    return meshAName + name1 + meshBName + name2;
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

void Foam::regionInterface::updateInterpolatorAndGlobalPatches()
{
    Info << "Updating interpolator and global patches" << endl;

    if (interfaceToInterfacePtr_.empty())
    {
        interfaceToInterface();
    }
    else if (interpolatorUpdateFrequency_ != 0)
    {
        if
        (
            ((runTime().timeIndex() - 1) % interpolatorUpdateFrequency_) == 0
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
Foam::regionInterface::interfaceToInterface() const
{
    if (interfaceToInterfacePtr_.empty())
    {
        makeInterfaceToInterface();
    }

    return interfaceToInterfacePtr_();
}

void Foam::regionInterface::attach()
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

void Foam::regionInterface::detach()
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


void Foam::regionInterface::updateUs()
{
    Us().internalField() = Up();

    correctUsBoundaryConditions();
}


void Foam::regionInterface::updatePhis()
{
    Phis() = fac::interpolate(Us()) & aMesh().Le();
}


void Foam::regionInterface::updateK()
{
    areaScalarField& curv = 
        const_cast<areaScalarField&>
        (
           aMesh().faceCurvatures()
        );

    correctCurvature(curv);

    curv.correctBoundaryConditions();

    K() == aMesh().faceCurvatures();


//    K().internalField() = 
//        const_cast<areaScalarField&>
//        (
//           aMesh().faceCurvatures()
//        );

//    correctCurvature(K());

//    K().correctBoundaryConditions();

//    if (runTime().outputTime())
//    {
//        K().write();
//    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "regionInterfaceTemplates.C"
#endif


// ************************************************************************* //
