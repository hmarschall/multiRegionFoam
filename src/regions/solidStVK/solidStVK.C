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
#include "transformGeometricField.H"


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

// * * * * * * * * * * * * Privat Member Functions * * * * * * * * * * * * * //

void  Foam::regionTypes::solidStVK::correctSigma(volSymmTensorField& sigma)
{
    // Calculate the right Cauchy–Green deformation tensor
    const volSymmTensorField c(symm(F_().T() & F_()));

    // Calculate the Green strain tensor
    const volSymmTensorField E(0.5*(c - I));

    // Calculate the 2nd Piola Kirchhoff stress
    const volSymmTensorField S(2.0*mu_*E + lambda_*tr(E)*I);

    // Convert the 2nd Piola Kirchhoff stress to the Cauchy stress
    // sigma = (1.0/J)*symm(F() & S & F().T());
    sigma = (1.0/J_())*transform(F_(), S);
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

    F_(nullptr),
    Finv_(nullptr),
    J_(nullptr),
    impK_(nullptr),
    impKf_(nullptr),
    rImpK_(nullptr),

    D_(nullptr),
    DD_(nullptr),
    pointD_(nullptr),
    pointDD_(nullptr),
    U_(nullptr),
    gradD_(nullptr),
    gradDD_(nullptr),
    sigma_(nullptr),

    stabilisationPtr_()
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

    F_ = lookupOrRead<volTensorField>
    (
        mesh(),
        "F",
        dimensionedTensor("I", dimless, I)
    );

    Finv_ = lookupOrRead<volTensorField>
    (
        mesh(),
        "Finv",
        dimensionedTensor("I", dimless, I)
        // false,
        // false,
        // inv(F_)
    );

    J_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "J",
        dimensionedScalar("1", dimless, 1.0)
        // false,
        // false,
        // det(F_)
    );

    impK_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "impK",
        (2.0*mu_ + lambda_)
    );

    impKf_ = lookupOrRead<surfaceScalarField>
    (
        mesh(),
        "impKf",
        (2.0*mu_ + lambda_)
    );

    rImpK_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "rImpK",
        false,
        false,
        1.0/impK_()
    );

    D_ = lookupOrRead<volVectorField>
    (
        mesh(),
        "D"
    );

    DD_ = lookupOrRead<volVectorField>
    (
        mesh(),
        "DD",
        dimensionedVector("0", dimLength, vector::zero)
    );

    gradD_ = lookupOrRead<volTensorField>
    (
        mesh(),
        "gradD",
        // true,
        // false,
        // fvc::grad(D_())
        dimensionedTensor("0", dimless, tensor::zero)
    );

    gradDD_ = lookupOrRead<volTensorField>
    (
        mesh(),
        "gradDD",
        // true,
        // false,
        // fvc::grad(DD_())
        dimensionedTensor("0", dimless, tensor::zero)
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
            dimensionedVector("0", vector::zero)
        )
    );
    cpi_.interpolate(D_(), pointD_());

    pointDD_.reset
    (
        new pointVectorField
        (
            IOobject
            (
                "pointDD",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            pointMesh::New(mesh()),
            dimensionedVector("0", vector::zero)
        )
    );
    cpi_.interpolate(DD_(), pointDD_());

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

    // Force old time fields to be stored
    D_().oldTime().oldTime();
    DD_().oldTime().oldTime();
    pointD_().oldTime();
    pointDD_().oldTime();
    gradD_().oldTime();
    gradDD_().oldTime();
    sigma_().oldTime();

    // Create stabilisation object
    if (!mechanicalProperties_.found("stabilisation"))
    {
        // If the stabilisation sub-dict is not found, we will add it with
        // default settings
        dictionary stabDict;
        stabDict.add("type", "RhieChow");
        stabDict.add("scaleFactor", 0.1);
        mechanicalProperties_.add("stabilisation", stabDict);
    }

    stabilisationPtr_.set
    (
        new momentumStabilisation
        (
            mechanicalProperties_.subDict("stabilisation")
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::solidStVK::~solidStVK()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::solidStVK::correct()
{
    // pointD_().storePrevIter();

    D_().correctBoundaryConditions();

    // DPrevOuterIter_() = D_();
}


Foam::scalar Foam::regionTypes::solidStVK::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::solidStVK::setCoupledEqns()
{
    D_().storePrevIter();

    DEqn =
    (
        rho_*fvm::d2dt2(D_())
     ==
        fvm::laplacian(impKf_(), D_(), "laplacian(DD,D)")
      - fvc::laplacian(impKf_(), D_(), "laplacian(DD,D)")
      + fvc::div(J_()*Finv_() & sigma_())
      + stabilisationPtr_().stabilisation(DD_(), gradDD_(), impK_)
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

    // Increment of displacement
    DD_() = D_() - D_().oldTime();

    // Update gradient of displacement
    gradD_() = fvc::grad(D_());

    // Update gradient of displacement increment
    gradDD_() = gradD_() - gradD_().oldTime();

    // Total deformation gradient
    F_() = I + gradD_().T();

    // Inverse of the deformation gradient
    Finv_() = inv(F_());

    // Jacobian of the deformation gradient
    J_() = det(F_());

    // Calculate the cauchy stress
    correctSigma(sigma_());

    // Interpolate cell displacements to vertices
    cpi_.interpolate(D_(), pointD_());

    // Increment of point displacement
    pointDD_() = pointD_() - pointD_().oldTime();

    // Velocity
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
