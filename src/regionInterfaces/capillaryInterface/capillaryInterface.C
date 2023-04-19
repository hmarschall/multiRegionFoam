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

#include "capillaryInterface.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionInterfaces
{
    defineTypeNameAndDebug(capillaryInterface, 0);

    addToRunTimeSelectionTable
    (
        regionInterface,
        capillaryInterface,
        IOdictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionInterfaces::capillaryInterface::capillaryInterface
(
    const word& type,
    const dictionary& dict,
    const Time& runTime,
    const fvPatch& patchA,
    const fvPatch& patchB
)
:
    regionInterface(type, dict, runTime, patchA, patchB),

    dict_(dict),

    sigma0_
    (
        dict_.lookup("sigma")
    ),
    sigma_
    (
        areaScalarField
        (
            IOobject
            (
                "sigma" + aMesh().mesh().name() + patchA.name(),
                runTime.timeName(),
                aMesh().thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh(),
            sigma0_,
            zeroGradientFaPatchScalarField::typeName
        )
    ),

    UsPtr_(),
    phisPtr_()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionInterfaces::capillaryInterface::clearOut() const
{
    UsPtr_.clear();
    phisPtr_.clear();

    regionInterface::clearOut();
}

void Foam::regionInterfaces::capillaryInterface::makeUs() const
{
    if (!UsPtr_.empty())
    {
        FatalErrorIn("regionInterface::makeUs()")
            << "surface velocity field already exists"
            << abort(FatalError);
    }

    // Set patch field types for Us
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

    // Set surface velocity
    UsPtr_.reset
    (
        new areaVectorField
        (
            IOobject
            (
                patchA().name() + "Us",
                runTime().timeName(),
                meshA(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh(),
            dimensioned<vector>("Us", dimVelocity, vector::zero),
            patchFieldTypes
        )
    );
}

void Foam::regionInterfaces::capillaryInterface::makePhis() const
{
    if (!phisPtr_.empty())
    {
        FatalErrorIn("regionInterface::makePhis()")
            << "surface fluid flux already exists"
            << abort(FatalError);
    }

    phisPtr_.reset
    (
        new edgeScalarField
        (
            IOobject
            (
                patchA().name() + "Phis",
                runTime().timeName(),
                meshA(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            linearEdgeInterpolate(Us()) & aMesh().Le()
        )
    );
}

void Foam::regionInterfaces::capillaryInterface::correctUsBoundaryConditions()
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

void Foam::regionInterfaces::capillaryInterface::updateUs()
{
    if (!meshA().foundObject<volVectorField>("U"))
    {
        return;
    }

    const volVectorField& U = meshA().lookupObject<volVectorField>("U");

    const fvBoundaryMesh& fvbm = meshA().boundary();

    const fvPatch& p = fvbm[patchAID()];

    Us().internalField() = p.lookupPatchField<volVectorField, vector>(U.name());

    correctUsBoundaryConditions();
}

void Foam::regionInterfaces::capillaryInterface::updatePhis()
{
    Phis() = fac::interpolate(Us()) & aMesh().Le();
}

tmp<vectorField>
Foam::regionInterfaces::capillaryInterface::surfaceTensionForce() const
{
    return
    (
        sigma0_.value()
       *fac::edgeIntegrate
        (
            aMesh().Le()
            *aMesh().edgeLengthCorrection()
        )().internalField()
    );
}

tmp<vectorField>
Foam::regionInterfaces::capillaryInterface::tangentialSurfaceTensionForce() const
{
    vectorField nA = meshA().boundary()[patchAID()].nf();

    return
    (
        surfaceTensionForce()
      - normalSurfaceTensionForce()
    );
}

tmp<vectorField>
Foam::regionInterfaces::capillaryInterface::normalSurfaceTensionForce() const
{
    vectorField nA = meshA().boundary()[patchAID()].nf();

    return
    (
        sigma_*aMesh().faceCurvatures().internalField()*nA
    );
}

void Foam::regionInterfaces::capillaryInterface::correct()
{
    // Update transport properties
    updateUs();
    updatePhis();

    // Update interface physics
    // TODO: call function to calculate new sigma_ with
    // interface equation of state for contaminated surfaces
}

Foam::scalar Foam::regionInterfaces::capillaryInterface::getMinDeltaT()
{
    scalar minDeltaT = GREAT;

    if (meshA().foundObject<volScalarField>("rho"))
    {
        const scalarField& dE = aMesh().lPN();
        scalar minDE = gMin(dE);
        const scalarField& rhoA = meshA().lookupObject<volScalarField>("rho");
        scalar minRhoA = gMin(rhoA);

        scalar maxCapillaryCo =
            runTime().controlDict().lookupOrDefault<scalar>("maxCapillaryCo", 1.0);

        minDeltaT =
            maxCapillaryCo*
            sqrt
            (
                minRhoA*minDE*minDE*minDE/
                2.0/M_PI/(sigma0_.value() + SMALL)
            );

    }

    return minDeltaT;
}

Foam::vector Foam::regionInterfaces::capillaryInterface::totalSurfaceTensionForce() const
{
    const scalarField& S = aMesh().S();

    return gSum
        (
            S
           *sigma0_.value()
           *fac::edgeIntegrate
            (
                aMesh().Le()
                *aMesh().edgeLengthCorrection()
            )().internalField()
        );
}

Foam::vector Foam::regionInterfaces::capillaryInterface::totalViscousForce() const
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    vectorField nGradU =
        meshA().lookupObject<volVectorField>("U")
        .boundaryField()[patchAID()]
        .snGrad();

    dimensionedScalar muA
        (
            meshA().lookupObject<IOdictionary>("transportProperties")
            .lookup("mu")
        );

    vectorField viscousForces =
      - muA.value()*S
       *(
            nGradU
          + (fac::grad(Us())().internalField()&n)
          - (n*fac::div(Us())().internalField())
        );

    return gSum(viscousForces);
}

Foam::vector Foam::regionInterfaces::capillaryInterface::totalPressureForce() const
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    const scalarField& p =
        meshA().lookupObject<volScalarField>("p").boundaryField()[patchAID()];

    vectorField pressureForces = S*p*n;

    return gSum(pressureForces);
}



void Foam::regionInterfaces::capillaryInterface::surfaceForces() const
{
    Info << "Total surface tension force: "
         << totalSurfaceTensionForce() << endl;

    vector totalForce = totalViscousForce() + totalPressureForce();

    Info << "Total force: " << totalForce << endl;
}

void Foam::regionInterfaces::capillaryInterface::maxCourantNumber() const
{
    const scalarField& dE = aMesh().lPN();

    dimensionedScalar rhoA
        (
            meshA().lookupObject<IOdictionary>("transportProperties")
            .lookup("rho")
        );

    scalar CoNum = gMax
    (
        runTime().deltaT().value()/
        sqrt
        (
            rhoA.value()*dE*dE*dE/
            2.0/M_PI/(sigma0_.value() + SMALL)
        )
    );

    Info << "Max surface Courant Number = " << CoNum << endl;
}

void Foam::regionInterfaces::capillaryInterface::curvature() const
{
    const scalarField& K = aMesh().faceCurvatures().internalField();

    Info << "Free surface curvature: min = " << gMin(K)
        << ", max = " << gMax(K)
        << ", average = " << gAverage(K) << endl << flush;
}

void Foam::regionInterfaces::capillaryInterface::info() const
{
    // Curvature
    curvature();

    // Surface courant number
    maxCourantNumber();

    // Surface Forces
    surfaceForces();
}

// ************************************************************************* //
