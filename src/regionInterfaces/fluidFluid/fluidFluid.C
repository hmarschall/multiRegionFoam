/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "fluidFluid.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

#include "wedgeFaPatch.H"
#include "wallFvPatch.H"
#include "wedgeFaPatchFields.H"
#include "slipFaPatchFields.H"
#include "fixedValueFaPatchFields.H"
#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionInterfaces
{
    defineTypeNameAndDebug(fluidFluid, 0);

    addToRunTimeSelectionTable
    (
        regionInterface, 
        fluidFluid, 
        IOdictionary 
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
    
Foam::regionInterfaces::fluidFluid::fluidFluid
( 
    const Time& runTime,   
    const fvPatch& patchA, 
    const fvPatch& patchB  
)
:
    regionInterface(runTime, patchA, patchB),
    
    //runTime_(runTime),
    
    //patchA_(patchA),
    //patchB_(patchB),
    
    U_
    (
        meshA().lookupObject<volVectorField>("U") 
    ),    
    phi_
    (
        meshA().lookupObject<surfaceScalarField>("phi")
    ),    
    rhoA_
    (
        this->subDict(meshA().name()).lookup("rho")
    ),
    rhoB_
    (
        this->subDict(meshB().name()).lookup("rho")
    ),
    muA_
    (
        this->subDict(meshA().name()).lookup("mu")
    ),
    muB_
    (
        this->subDict(meshB().name()).lookup("mu")
    ),
    sigma0_
    (
        this->lookup("cleanSurfaceTension")
    ),
    g_
    (
        this->lookup("g")
    ),
    
    UsPtr_(NULL), 
    KPtr_(NULL),
    phisPtr_(NULL), 
    
    massFluxCorr_
    (
        meshA().boundaryMesh()[patchAID()].size(), 0
    ),
    correctContinuityError_
    (
        Switch
        (
            this->lookupOrDefault<Switch>("correctContinuityError", true)
        )
    ),
    curTimeIndex_(0),
    initialTotalVolume_(0.),
    totalVolume_(0.)

{
    // add
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::regionInterfaces::fluidFluid::clearOut()
{
    deleteDemandDrivenData(UsPtr_);
    deleteDemandDrivenData(KPtr_);
    deleteDemandDrivenData(phisPtr_);
}

void Foam::regionInterfaces::fluidFluid::makeK() const
{
    if (KPtr_)
    {
        FatalErrorIn("fluidFluid::makeK()")
            << "surface curvature field already exists"
            << abort(FatalError);
    }

    KPtr_ = new areaScalarField
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
    );
}


void Foam::regionInterfaces::fluidFluid::makePhis() const
{
    if (phisPtr_)
    {
        FatalErrorIn("fluidFluid::makePhis()")
            << "surface fluid flux already exists"
            << abort(FatalError);
    }

    phisPtr_ = new edgeScalarField
    (
        IOobject
        (
            phi_.name() + "s",
            runTime().constant(), 
            runTime(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearEdgeInterpolate(Us()) & aMesh().Le()
    );
} 


void Foam::regionInterfaces::fluidFluid::makeUs() const
{
    if (UsPtr_)
    {
        FatalErrorIn("fluidFluid::makeUs()")
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
    // ... lengthy code commented
    
    UsPtr_ = new areaVectorField
    (
        IOobject
        (
            U_.name() + "s",
            runTime().constant(), 
            runTime(), 
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh(),
        dimensioned<vector>("Us", dimVelocity, vector::zero),
        patchFieldTypes
    );   
    // ... lengthy code commented
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionInterfaces::fluidFluid::~fluidFluid()
{
    clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const vectorField& Foam::regionInterfaces::fluidFluid::Up()
{
    const fvBoundaryMesh& fvbm = meshA().boundary(); 

    const fvPatch& p = fvbm[patchAID()];

    return p.lookupPatchField<volVectorField, vector>(U_.name());
} 


void Foam::regionInterfaces::fluidFluid::correctUsBoundaryConditions()
{   
    forAll(Us().boundaryField(), patchI)
    {
        if
        (
            (*UsPtr_).boundaryField()[patchI].type()
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
                        U_.boundaryField()[ngbPolyPatchID].type()
                     == slipFvPatchVectorField::typeName
                    )
                 ||
                    (
                        U_.boundaryField()[ngbPolyPatchID].type()
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


void Foam::regionInterfaces::fluidFluid::updateUs()
{
    Us().internalField() = Up();

    correctUsBoundaryConditions();
}


void Foam::regionInterfaces::fluidFluid::updateSurfaceFlux()
{
    Phis() = fac::interpolate(Us()) & aMesh().Le();
}


void Foam::regionInterfaces::fluidFluid::calcCurvatureAxis
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
    // ... lengthy code commented
    }
}


void Foam::regionInterfaces::fluidFluid::updateK()
{
    K().internalField() = 
        const_cast<areaScalarField&>
        (
           aMesh().faceCurvatures()
        );

    calcCurvatureAxis(K());

    K().correctBoundaryConditions();
}


tmp<scalarField> Foam::regionInterfaces::fluidFluid::massFluxCorr() const
{
    tmp<scalarField> pf(new scalarField(massFluxCorr_));

    if (!correctContinuityError_ || curTimeIndex_ == runTime().timeIndex())
    {
        return pf;
    }
    else
    {
        initialTotalVolume_ = sum(meshA().cellVolumes()); 
        totalVolume_ = sum(meshA().cellVolumes()); 

        curTimeIndex_ = runTime().timeIndex();
    }

    scalar newTotalVolume = sum(meshA().cellVolumes()); 

    Info<< "Volume: new = " << newTotalVolume << " old = " << totalVolume_
        << " relative change = " << Foam::mag(newTotalVolume - totalVolume_) 
        << " global error = " << scalar(1-newTotalVolume/initialTotalVolume_) 
        << endl;

    //TODO: does not work, needs further dev
    totalVolume_ = newTotalVolume;
    
    // why redefining phi, U, rho
    const surfaceScalarField& phi = 
        meshA().objectRegistry::lookupObject<surfaceScalarField>("phi");
        
    const volVectorField& U = 
        meshA().objectRegistry::lookupObject<volVectorField>("U");

    const volScalarField& rho = 
        meshA().objectRegistry::lookupObject<volScalarField>("rho");

    const fvPatchScalarField& rhop = 
        rho.boundaryField()[patchAID()];

    scalar fluxCorr = 0.;

    scalar massIn = 0.;
    scalar massOut = 0.;

    forAll (phi.boundaryField(), patchi)
    {
        const fvPatchVectorField& Up = U.boundaryField()[patchi];

        const fvPatchScalarField& rhop = 
            rho.boundaryField()[patchi];

        fvsPatchScalarField& phip = const_cast<fvsPatchScalarField&>
        (
            phi.boundaryField()[patchi]
        );

        if 
        (
            !Up.coupled()
        )
        {
            forAll (phip, i)
            {
                if (phip[i] < 0.0)
                {
                    massIn -= rhoA_.value()*phip[i];
                }
                else
                {
                    massOut += rhoA_.value()*phip[i];
                }
            }
        }
    }

    Info << "  massIn - massOut : " << (massIn - massOut) << endl;

    scalar dt = runTime().deltaT().value(); 

    scalar sumV = 0.;
    scalar sumV0 = 0.;
    scalar sumV00 = 0.;

    if (runTime().timeIndex() > 3) 
    {
        scalarField V = meshA().V().field(); 
        scalarField V0 = meshA().V0().field();
        scalarField V00 = meshA().V00().field();

        sumV = gSum(V);
        sumV0 = gSum(V0);
        sumV00 = gSum(V00);

        Info<< "  Volume continuity errors : "
            << "volume = " << sumV
            << ", old volume = " << sumV0
            << ", global = " << (sumV - sumV0) << endl;
    }

    reduce(sumV, sumOp<scalar>());
    reduce(sumV0, sumOp<scalar>());
    reduce(sumV00, sumOp<scalar>());


    scalarField magSf = meshA().boundary()[patchAID()].magSf(); 
    scalar magSftot = gSum(magSf);
    reduce(magSftot, sumOp<scalar>());

    Info << "  Interfacial area : " << magSftot << endl;

    scalar relax = this->lookupOrDefault<scalar>("relax", 0.1);

    pf() = rhoB_.value()*(sumV - sumV0)/(dt*magSftot);

    Info << "  Avrg mass flux correction : " << average(pf()) << endl;

    return pf;
}


// virtual functions from regionInterface.H 
// but already defined in regionIinterface.C

void Foam::regionInterfaces::fluidFluid::attach() 
{
    // do nothing, add as required
}

void Foam::regionInterfaces::fluidFluid::detach() 
{
    // do nothing, add as required
}

void Foam::regionInterfaces::fluidFluid::updateInterpolatorAndGlobalPatches() 
{
    // do nothing, add as required
}

// ************************************************************************* //
