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
#include "ionomer.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "ionomerParameters.H"
#include "globalFCParameters.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(ionomer, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        ionomer,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::regionTypes::ionomer::updateIonomerProperties()
{
    // Diffusion coefficient of water through ionomer
    DLambda_() = (3.842*pow(lambda_(),3) - 32.03*sqr(lambda_()) + 67.75*lambda_())/(pow(lambda_(),3) - 2.115*sqr(lambda_()) - 33.013*lambda_() + 103.37)*DLambda0_*exp(ELambda_/RGas_*((1/TRef_) - (1/T_())));

    // electro-osmotic drag coefficient
    xi_ = 2.5*lambda_()/22;

    // volume fraction water
    f_ = lambda_()*VW_/(lambda_()*VW_ + VM_);
    fCond_ = fComp_;

    // protonic conductivity
    kappa_() = pow(pos(fCond_)*fCond_,1.5)*kappa0_*exp(EKappa_/RGas_*((1/TRef_) - (1/T_())));
}

void Foam::regionTypes::ionomer::updateSourceTerms()
{
    // heat source - joule heating protons
    //sT_ = kappa_()*(fvc::grad(phiP_())&fvc::grad(phiP_()));

    // water content source - electro-osmotic drag
    sLambda_ = fvc::laplacian(xi_*kappa_()/FConst_, phiP_());
}	

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::ionomer::ionomer
(
    const Time& runTime,
    const word& regionName
)
:
    regionType(runTime, regionName),

    regionName_(regionName),

    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            mesh().time().constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    materialProperties_
    (
        IOobject
        (
            "materialProperties",
            mesh().time().constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    operatingConditions_
    (
        IOobject
        (
            "operatingConditions",
            mesh().time().constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    cv_(transportProperties_.lookup("cv")),
    rho_(transportProperties_.lookup("rho")),
    k_(nullptr),
    kappa_(nullptr),
    VW_(materialProperties_.lookup("VW")),
    VM_(materialProperties_.lookup("VM")),
    RH_(operatingConditions_.lookup("RH")),
    DLambda_(nullptr),
    f_
    (
        IOobject
        (
           "f",
           mesh().time().timeName(),
           mesh(),
           IOobject::READ_IF_PRESENT,
           IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("f0", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.26)
    ),
    fCond_
    (
        IOobject
        (
            "fCond",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("f0_", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.2)
    ),
    xi_
    (
        IOobject
        (
            "xi",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("xi0", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.1)
    ),
    sT_
    (
        IOobject
        (
            "sT",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("sT0", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0)
    ),
    sLambda_
    (
        IOobject
        (
            "sLambda",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("sLambda0", dimensionSet(0, -3, -1, 0, 1, 0, 0), 0)
    ),
    T_(nullptr),
    phiP_(nullptr),
    lambda_(nullptr)
{
    k_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "k",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar(transportProperties_.lookup("k"))
        )
    );

    kappa_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "kappa",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("kappaInit", dimensionSet(-1, -3, 3, 0, 0, 2, 0), 3)
        )
    );

    DLambda_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "DLambda",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            DLambda0_
        )
    );

    T_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "T",
                mesh().time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );

    phiP_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "phiP",
                mesh().time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );

    lambda_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "lambda",
                mesh().time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );

    // thermal conductivity
    k_() = dimensionedScalar(transportProperties_.lookup("k"));

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::ionomer::~ionomer()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::ionomer::correct()
{
    // do nothing, add as required
}


void Foam::regionTypes::ionomer::setRDeltaT()
{
    // do nothing, add as required
}


void Foam::regionTypes::ionomer::setCoupledEqns()
{
    if (runTime().timeIndex() != 0)
    {
        // update variables
        // ionomer properties
        updateIonomerProperties();

        // source terms
        updateSourceTerms();
    }

    // fourier heat conduction
    fvScalarMatrix TEqn =
    (
        rho_*cv_*fvm::ddt(T_())
     ==
        fvm::laplacian(k_(), T_(), "laplacian(k,T)")
       //+sT_
    );

    // ohm's law for protons
    fvScalarMatrix phiPEqn =
    (
        -fvm::laplacian(kappa_(), phiP_(), "laplacian(kappa,phiP)")
    );

    // water transport in ionomer
    fvScalarMatrix lambdaEqn =
    (
        1/VM_*fvm::ddt(lambda_())
       ==
        fvm::laplacian(DLambda_()/VM_, lambda_(), "laplacian(DLambda,lambda)")
        +sLambda_
    );
    
    fvScalarMatrices.set
    (
        T_().name() + mesh().name() + "Eqn",
        new fvScalarMatrix(TEqn)
    );

    fvScalarMatrices.set
    (
        phiP_().name() + mesh().name() + "Eqn",
        new fvScalarMatrix(phiPEqn)
    );

    fvScalarMatrices.set
    (
        lambda_().name() + mesh().name() + "Eqn",
        new fvScalarMatrix(lambdaEqn)
    );
}

void Foam::regionTypes::ionomer::updateFields()
{
    Info<< "Temperature = "
            << T_().weightedAverage(mesh().V()).value()
            << " Min(T) = " << min(T_()).value()
            << " Max(T) = " << max(T_()).value()
            << endl;

    Info<< "Water content = "
            << lambda_().weightedAverage(mesh().V()).value()
            << " Min(lambda) = " << min(lambda_()).value()
            << " Max(lambda) = " << max(lambda_()).value()
            << endl;

    Info<< "Electrolye Potential = "
            << phiP_().weightedAverage(mesh().V()).value()
            << " Min(phiP) = " << min(phiP_()).value()
            << " Max(phiP) = " << max(phiP_()).value()
            << endl;
}

void Foam::regionTypes::ionomer::solveRegion()
{
    // do nothing, add as required
}


// ************************************************************************* //
