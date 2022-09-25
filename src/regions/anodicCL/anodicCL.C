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
#include "anodicCL.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "anodicCLParameters.H"
#include "globalFCParameters.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(anodicCL, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        anodicCL,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::regionTypes::anodicCL::updateIonomerProperties()
{
    // water diffusion through ionomer
    DLambda_() = pow(epsilonI_,1.5)*(3.842*pow(lambda_(),3) - 32.03*sqr(lambda_()) + 67.75*lambda_())/(pow(lambda_(),3) - 2.115*sqr(lambda_()) - 33.013*lambda_() + 103.37)*DLambda0_*exp(ELambda_/RGas_*((1/TRef_) - (1/T_())));

    // electro-osmotic drag coefficient
    xi_ = 2.5*lambda_()/22;

    // volume fraction water
    f_ = lambda_()*VW_/(lambda_()*VW_ + VM_);
    fCond_ = f_ - fComp_;

    // protonic conductivity
    kappa_() = pow(pos(fCond_)*fCond_,1.5)*kappa0_*pow(epsilonI_,1.5)*exp(EKappa_/RGas_*((1/TRef_) - (1/T_())));
}

void Foam::regionTypes::anodicCL::updateGasSpeciesTransportProperties()
{
    // ideal gas concentration
    c_ = p_/(RGas_*T_());

    // effective diffusion coefficient hydrogen
    DEffH2_() = epsilonP_/sqr(tau_)*DH2_*pow((T_()/TRef_),1.5)*(pRef_/p_);

    // effective diffusion coefficient vapor
    DEffV_() = epsilonP_/sqr(tau_)*DV_*pow((T_()/TRef_),1.5)*(pRef_/p_);
}

void Foam::regionTypes::anodicCL::updateAbsorptionDesorption()
{
    // saturation vapor fraction
    xVSat_ = (exp(23.1963 - (TRefP1_/(T_() - TRefP2_)))*pDim_)/p_;

    // equilibrium water content in ionomer
    lambdaEq_ = 0.043 + (17.81*xV_()/xVSat_) - (39.85*pow((xV_()/xVSat_),2)) + (36*pow((xV_()/xVSat_),3));

    // ab-/desorption rate
    kSorp_ = (pos(lambda_() - lambdaEq_)*aD_ + (1 - pos(lambda_() - lambdaEq_))*aA_)*f_*exp((ELambda_/RGas_)*((1/TRef_) - (1/T_())));
}

void Foam::regionTypes::anodicCL::updateElectrochemistry()
{
    // reaction exchange current density
    j0_ = jStern_*exp(ER_/RGas_*(1/TRef_ - 1/T_()));

    // galvani potential difference
    deltaPhi_ = phiE_() - phiP_();

    // reversible potential difference
    deltaPhi0_ = -(T_()*deltaS_/(2*FConst_)) - (RGas_*T_()/(2*FConst_))*log(xH2_()*p_/pRef_);

    // overpotential
    eta_ = deltaPhi_ - deltaPhi0_;

    // volumetric exchange current density (Butler-Volmer)
    j_ = -j0_*a_*(exp(2*beta_*FConst_*eta_/RGas_/T_()) - exp(-2*(1 - beta_)*FConst_*eta_/RGas_/T_()));
}

void Foam::regionTypes::anodicCL::updateSourceTerms()
{
    // heat Source - joule heating electrons & protons, sorption heat, reaction heat
    sT_ = (sigma_()*(fvc::grad(phiE_())&fvc::grad(phiE_())) + kappa_()*(fvc::grad(phiP_())&fvc::grad(phiP_())) + (kSorp_/d_/VM_)*(lambdaEq_ - lambda_())*HAD_ + j_*eta_ - (j_/2/FConst_)*T_()*deltaS_)/T_();

    // mass source water content in ionomer
    sLambda_ = ((kSorp_/d_/VM_)*(lambdaEq_ - lambda_()))/lambda_();

    // mass source vapor / phase change & sorption
    sV_ = - ((kSorp_/d_/VM_)*(lambdaEq_ - lambda_()))/xV_();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::anodicCL::anodicCL
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
    electrochemicalProperties_
    (
        IOobject
        (
            "electrochemicalProperties",
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
    kappa_(nullptr),
    DH2_(transportProperties_.lookup("DH2")),
    DV_(transportProperties_.lookup("DV")),
    epsilonP_(materialProperties_.lookup("epsilonP")),
    tau_(materialProperties_.lookup("tau")),
    epsilonI_(materialProperties_.lookup("epsilonI")),
    d_(materialProperties_.lookup("d")),
    MW_(materialProperties_.lookup("MW")),
    VW_(materialProperties_.lookup("VW")),
    VM_(materialProperties_.lookup("VM")),
    a_(electrochemicalProperties_.lookup("a")),
    beta_(electrochemicalProperties_.lookup("beta")),
    deltaS_(electrochemicalProperties_.lookup("deltaS")),
    ER_(electrochemicalProperties_.lookup("ER")),
    jStern_(electrochemicalProperties_.lookup("jStern")),
    p_(operatingConditions_.lookup("p")),
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
        dimensionedScalar("f0", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.2)
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
    lambdaEq_
    (
        IOobject
        (
            "lambdaEq",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("lambdaEq0", dimensionSet(0, 0, 0, 0, 0, 0, 0), 12)
    ),
    kSorp_
    (
        IOobject
        (
            "kSorp",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("kSorp0", dimensionSet(0, 1, -1, 0, 0, 0, 0), 0)
    ),
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
    xVSat_
    (
        IOobject
        (
            "xVSat",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("xVSat0", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.2)
    ),
    j0_
    (
        IOobject
        (
            "j0",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        jStern_
    ),
    deltaPhi_
    (
        IOobject
        (
            "deltaPhi",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("deltaPhiInit", dimensionSet(1, 2, -3, 0, 0, -1, 0), 1)
    ),
    deltaPhi0_
    (
        IOobject
        (
            "deltaPhi0",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("deltaPhi0Init", dimensionSet(1, 2, -3, 0, 0, -1, 0), 1.19)
    ),
    eta_
    (
        IOobject
        (
            "eta",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("eta0", dimensionSet(1, 2, -3, 0, 0, -1, 0), 0.1)
    ),
    j_
    (
        IOobject
        (
            "j",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("jInit", dimensionSet(0, -3, 0, 0, 0, 1, 0), 0)
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
        dimensionedScalar("sT0", dimensionSet(1, -1, -3, -1, 0, 0, 0), 0)
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
    sV_
    (
        IOobject
        (
            "sV",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("sV0", dimensionSet(0, -3, -1, 0, 1, 0, 0), 0)
    ),
    T_(nullptr),
    phiE_(nullptr),
    phiP_(nullptr),
    lambda_(nullptr),
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
    // electric conductivity
    sigma_() = dimensionedScalar(transportProperties_.lookup("sigma"));

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::anodicCL::~anodicCL()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::anodicCL::correct()
{
    // do nothing, add as required
}


void Foam::regionTypes::anodicCL::setRDeltaT()
{
    // do nothing, add as required
}


void Foam::regionTypes::anodicCL::setCoupledEqns()
{
    
    if (runTime().timeIndex() != 0)
    {
        // update fields
        // ionomer properties
        updateIonomerProperties();

        // gas species transport
        updateGasSpeciesTransportProperties();

        // ab-/desorption
        updateAbsorptionDesorption();

        // electrochemistry
        updateElectrochemistry();

        // source terms
        updateSourceTerms();
    }

    // set anodic Eqns
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
        - fvm::laplacian(sigma_(), phiE_(), "laplacian(sigma,phiE)")
        ==
        - fvm::SuSp(-j_/phiE_(), phiE_())
    );

    // ohm's law for protons
    fvScalarMatrix phiPEqn =
    (
          fvm::laplacian(kappa_(), phiP_(), "laplacian(kappa,phiP)")
        ==
        - fvm::SuSp(-j_/phiP_(), phiP_())
    );

    // water transport in ionomer
    fvScalarMatrix lambdaEqn =
    (
          1/VM_*fvm::ddt(lambda_())
        - fvm::laplacian(DLambda_()/VM_, lambda_(), "laplacian(DLambda,lambda)")
        + fvc::laplacian(xi_*kappa_()/FConst_, phiP_(), "laplacian(kappa,phiP)")
        ==
          fvm::SuSp(sLambda_, lambda_()) 
    );

    // fick diffusion for vapor
    fvScalarMatrix xVEqn =
    (
          c_*fvm::ddt(xV_())
        - c_*fvm::laplacian(DEffV_(), xV_(), "laplacian(D,x)")
        ==
          fvm::SuSp(sV_, xV_())
    );

    // fick diffusion for hydrogen
    fvScalarMatrix xH2Eqn =
    (
          c_*fvm::ddt(xH2_())
        - c_*fvm::laplacian(DEffH2_(), xH2_(), "laplacian(D,x)")
        ==
          fvm::SuSp(-j_/2/FConst_/xH2_(), xH2_())
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
        phiP_().name() + mesh().name() + "Eqn",
        new fvScalarMatrix(phiPEqn)
    );

    fvScalarMatrices.set
    (
        lambda_().name() + mesh().name() + "Eqn",
        new fvScalarMatrix(lambdaEqn)
    );

    fvScalarMatrices.set
    (
        xV_().name() + mesh().name() + "Eqn",
        new fvScalarMatrix(xVEqn)
    );  
      
    fvScalarMatrices.set
    (
        xH2_().name() + mesh().name() + "Eqn",
        new fvScalarMatrix(xH2Eqn)
    );
}

void Foam::regionTypes::anodicCL::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::anodicCL::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::anodicCL::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::anodicCL::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::anodicCL::postSolve()
{
    // do nothing, add as required
}



// ************************************************************************* //
