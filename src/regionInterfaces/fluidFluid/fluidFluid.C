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

    transportPropertiesA_
    (
        IOobject
        (
            "transportProperties",
            fileName(runTime.caseConstant()/meshA().name()),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    transportPropertiesB_
    (
        IOobject
        (
            "transportProperties",
            fileName(runTime.caseConstant()/meshB().name()),
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
        transportPropertiesA_.lookup("rho")
    ),
    rhoB_
    (
        transportPropertiesB_.lookup("rho")
    ),
    muA_
    (
        transportPropertiesA_.lookup("mu")
    ),
    muB_
    (
        transportPropertiesB_.lookup("mu")
    ),
    sigma0_
    (
        interfaceProperties().subDict(name()).lookup("sigma")
    ),
    g_
    (
        gravitationalProperties_.lookup("g")
    )
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::regionInterfaces::fluidFluid::makeK() const
{
    if (!KPtr_.empty())
    {
        FatalErrorIn("fluidFluid::makeK()")
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


void Foam::regionInterfaces::fluidFluid::makePhis() const
{
    if (!phisPtr_.empty())
    {
        FatalErrorIn("fluidFluid::makePhis()")
            << "surface fluid flux already exists"
            << abort(FatalError);
    }

    phisPtr_.set
    (
        new edgeScalarField
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
        )
    );
} 


void Foam::regionInterfaces::fluidFluid::makeUs() const
{
    if (!UsPtr_.empty())
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
    
    UsPtr_.set
    (
        new areaVectorField
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
        )
    );
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

// ************************************************************************* //
