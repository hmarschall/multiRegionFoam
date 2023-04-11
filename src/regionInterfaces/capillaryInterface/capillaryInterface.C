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

    // Get data needed for Us calculation
    vectorField nA = meshA().boundary()[patchAID()].nf();

    const volVectorField& UA = meshA().lookupObject<volVectorField>("U");
    vectorField UAp = UA.boundaryField()[patchAID()].patchInternalField();

    const volVectorField& UB = meshB().lookupObject<volVectorField>("U");
    vectorField UBp = interpolateFacesFromB
        (
            UB.boundaryField()[patchBID()].patchInternalField()()
        );

    const surfaceScalarField& phi = meshA()
        .lookupObject<surfaceScalarField>("phi");
    scalarField phiP = phi.boundaryField()[patchAID()];

    dimensionedScalar muA
        (
            meshA().lookupObject<IOdictionary>("transportProperties")
            .lookup("mu")
        );

    dimensionedScalar muB
        (
            meshB().lookupObject<IOdictionary>("transportProperties")
            .lookup("mu")
        );

    scalarField DnA = meshA().boundary()[patchAID()].deltaCoeffs();

    scalarField DnB = interpolateFacesFromB
        (
            meshB().boundary()[patchBID()].deltaCoeffs()
        );


    // Calculate Us
    vectorField UtPA = UAp;
    UtPA -= nA*(nA & UtPA);

    vectorField UtPB = UBp;
    UtPB -= nA*(nA & UtPB);

    vectorField UnFs =nA*phiP/meshA().boundary()[patchAID()].magSf();

    Us().internalField() += UnFs - nA*(nA&Us().internalField());
    correctUsBoundaryConditions();

    vectorField UtFs =
        muA.value()*DnA*UtPA
      + muB.value()*DnB*UtPB
      + (muB.value() - muA.value())
       *(fac::grad(Us())&aMesh().faceAreaNormals())().internalField()
      + tangentialSurfaceTensionForce();

    UtFs /= muA.value()*DnA + muB.value()*DnB;

    Us().internalField() = UnFs + UtFs;

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
        sigma_
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


// ************************************************************************* //
