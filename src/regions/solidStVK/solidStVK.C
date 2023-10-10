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

#include "fvCFD.H"
#include "solidStVK.H"
#include "addToRunTimeSelectionTable.H"
#include "dimensionedTypes.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(solidStVK, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        solidStVK,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::solidStVK::solidStVK
(
    const Time& runTime,
    const word& regionName
)
:
    regionType(runTime, regionName),

    regionName_(regionName),

    cpi_
    (
        volPointInterpolation::New(mesh())
    ),

    mechanicalProperties_
    (
        IOobject
        (
            "mechanicalProperties",
            mesh().time().constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    rho_(dimensionedScalar(mechanicalProperties_.lookup("rho"))),
    E_(dimensionedScalar(mechanicalProperties_.lookup("E"))),
    nu_(dimensionedScalar(mechanicalProperties_.lookup("nu"))),
    mu_(E_/(2.0*(1.0 + nu_))),
    lambda_(nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))),
    threeK_(E_/(1.0 - 2.0*nu_)),

    D_(nullptr),
    pointD_(nullptr),
    U_(nullptr),
    gradD_(nullptr),
    sigma_(nullptr),
    DPrevOuterIter_(nullptr)
{
    // Correct Lamé coefficients for plane stress
    Switch planeStress(mechanicalProperties_.lookup("planeStress"));

    if (planeStress)
    {
        Info<< "Plane Stress\n" << endl;

        //- change lambda and threeK for plane stress
        lambda_ = nu_*E_/((1.0 + nu_)*(1.0 - nu_));
        threeK_ = E_/(1.0 - nu_);
    }

    //D_ = lookupOrRead<volVectorField>(mesh(), "D");

    D_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "D",
                runTime.timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );

    gradD_ = lookupOrRead<volTensorField>
    (
        mesh(),
        "gradD",
        // false,
        // true,
        //fvc::grad(D_())
        dimensionedTensor("0", dimless, tensor::zero),
        true
    );

    // Change lookupOrRead functionality to also handle pointMesh fields
    pointD_.reset
    (
        new pointVectorField
        (
            IOobject
            (
                "pointD",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            pointMesh::New(mesh()),
            dimensionedVector("", vector::zero)
        )
    );
    cpi_.interpolate(D_(), pointD_());
    pointD_->storePrevIter();

    U_ = lookupOrRead<volVectorField>
    (
        mesh(),
        "U",
        // false,
        // true,
        //fvc::ddt(D_())
        dimensionedVector("0", dimLength/dimTime, vector::zero),
        true
    );

    sigma_ = lookupOrRead<volSymmTensorField>
    (
        mesh(),
        "sigma",
        // false,
        // true,
        // symm
        // (
        //     (mu_*(gradD_() + gradD_().T()))
        //     + ((lambda_*I)*tr(gradD_()))
        //     + (mu_*(gradD_() & gradD_().T()))
        //     + (0.5*lambda_*I*tr(gradD_() & gradD_().T()))
        // )
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero),
        true
    );

    DPrevOuterIter_ = lookupOrRead<volVectorField>
    (
        mesh(),
        "DPrevOuterIter",
        false,
        true,
        D_()
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::solidStVK::~solidStVK()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::solidStVK::correct()
{
    pointD_().storePrevIter();

    D_().correctBoundaryConditions();

    DPrevOuterIter_() = D_();
}


Foam::scalar Foam::regionTypes::solidStVK::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::solidStVK::setCoupledEqns()
{

    #   include "solidStVKDdtCoeffs.H"

    D_().storePrevIter();

    DEqn =
    (
        rho_*(Cn*fvm::ddt(D_()) - Co*U_().oldTime() + Coo*U_().oldTime().oldTime())
     ==
        fvm::laplacian(2*mu_ + lambda_, D_(), "laplacian(DD,D)")
      - fvc::laplacian(2*mu_ + lambda_, D_(), "laplacian(DD,D)")
      + fvc::div(sigma_() & (I + gradD_()))
    );

    fvVectorMatrices.set
    (
        D_().name()
      + mesh().name() + "Mesh"
      + solidStVK::typeName + "Type"
      + "Eqn",
        &DEqn()
    );
}

void Foam::regionTypes::solidStVK::postSolve()
{
    D_().relax();

    gradD_() = fvc::grad(D_());
    sigma_() = symm((mu_*(gradD_() + gradD_().T()))
           + ((lambda_*I)*tr(gradD_()))
           + (mu_*(gradD_() & gradD_().T()))
           + (0.5*lambda_*I*tr(gradD_() & gradD_().T())));

    cpi_.interpolate(D_(), pointD_());

    pointD_().correctBoundaryConditions();

    U_() = fvc::ddt(D_());
}

void Foam::regionTypes::solidStVK::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::solidStVK::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::solidStVK::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::solidStVK::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::solidStVK::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //
