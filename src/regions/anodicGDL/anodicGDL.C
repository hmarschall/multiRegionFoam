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
#include "anodicGDL.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "anodicGDLParameters.H"
#include "globalFCParameters.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(anodicGDL, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        anodicGDL,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::regionTypes::anodicGDL::updateGasSpeciesTransportProperties()
{
    // ideal gas concentration
    c_ = p_/(RGas_*T_());

    // effective diffusion coefficient hydrogen
    DEffH2_() = epsilonP_/sqr(tau_)*DH2_*pow((T_()/TRef_),1.5)*(pRef_/p_);

    // effective diffusion coefficient vapor
    DEffV_() = epsilonP_/sqr(tau_)*DV_*pow((T_()/TRef_),1.5)*(pRef_/p_);
}

void Foam::regionTypes::anodicGDL::updateSourceTerms()
{
    // heat Source - joule heating electrons
    sT_ = (sigma_()*(fvc::grad(phiE_())&fvc::grad(phiE_())))/T_();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::anodicGDL::anodicGDL
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
    sigma_(nullptr),
    DH2_(transportProperties_.lookup("DH2")),
    DV_(transportProperties_.lookup("DV")),
    epsilonP_(materialProperties_.lookup("epsilonP")),
    tau_(materialProperties_.lookup("tau")),
    p_(operatingConditions_.lookup("p")),
    c_
    (
        IOobject
        (
            "c",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("c0", dimensionSet(0, -3, 0, 0, 1, 0, 0), 52)
    ),
    DEffH2_(nullptr),
    DEffV_(nullptr),
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
        dimensionedScalar("sT0", dimensionSet(1, -1, -3, -1, 0, 0, 0), 0)
    ),
    T_(nullptr),
    phiE_(nullptr),
    xH2_(nullptr),
    xV_(nullptr)
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

    sigma_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "sigma",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar(transportProperties_.lookup("sigma"))
        )
    );

    DEffH2_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "DEffH2",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            DH2_
        )
    );

    DEffV_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "DEffV",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            DV_
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

    phiE_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "phiE",
                mesh().time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );

    xH2_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "xH2",
                mesh().time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );

    xV_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "xV",
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
    // electric conducivity
    sigma_() = dimensionedScalar(transportProperties_.lookup("sigma"));
    
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::anodicGDL::~anodicGDL()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::anodicGDL::correct()
{
    // do nothing, add as required
}


void Foam::regionTypes::anodicGDL::setRDeltaT()
{
    // do nothing, add as required
}


void Foam::regionTypes::anodicGDL::setCoupledEqns()
{
    if (runTime().timeIndex() != 0)
    {
        // update fields
        // gas species transport
        updateGasSpeciesTransportProperties();
    
        // source terms
        updateSourceTerms();
    }

    // set Eqns
    // fourier heat conduction
    fvScalarMatrix TEqn =
    (
          rho_*cv_*fvm::ddt(T_())
        - fvm::laplacian(k_(), T_(), "laplacian(k,T)")
        ==
        - fvm::SuSp(-sT_, T_())
    );

    // ohm's law for electrons
    fvScalarMatrix phiEEqn =
    (
        -fvm::laplacian(sigma_(), phiE_(), "laplacian(sigma,phiE)")
    );

    // fick diffusion for hydrogen
    fvScalarMatrix xH2Eqn =
    (
          c_*fvm::ddt(xH2_())
        ==
          fvm::laplacian(c_*DEffH2_(), xH2_(), "laplacian(D,x)")
    );

    // fick diffusion for vapor
    fvScalarMatrix xVEqn =
    (
          c_*fvm::ddt(xV_())
        ==
          fvm::laplacian(c_*DEffV_(), xV_(), "laplacian(D,x)")
    );

    fvScalarMatrices.set
    (
        T_().name() + mesh().name() + "Eqn",
        new fvScalarMatrix(TEqn)
    );

    fvScalarMatrices.set
    (
        phiE_().name() + mesh().name() + "Eqn",
        new fvScalarMatrix(phiEEqn)
    );

    fvScalarMatrices.set
    (
        xH2_().name() + mesh().name() + "Eqn",
        new fvScalarMatrix(xH2Eqn)
    );

    fvScalarMatrices.set
    (
        xV_().name() + mesh().name() + "Eqn",
        new fvScalarMatrix(xVEqn)
    );
}

void Foam::regionTypes::anodicGDL::updateFields()
{
    
}

void Foam::regionTypes::anodicGDL::solveRegion()
{
    // do nothing, add as required
}


// ************************************************************************* //
