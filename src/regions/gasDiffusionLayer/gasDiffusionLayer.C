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
#include "gasDiffusionLayer.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "gasDiffusionLayerParameters.H"
#include "globalFCParameters.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(gasDiffusionLayer, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        gasDiffusionLayer,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::gasDiffusionLayer::gasDiffusionLayer
(
    const fvMesh& mesh,
    const word& regionName
)
:
    regionType(mesh, regionName),

    mesh_(mesh),
    regionName_(regionName),

    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            this->time().constant(),
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    materialProperties_
    (
        IOobject
        (
            "materialProperties",
            this->time().constant(),
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    operatingConditions_
    (
        IOobject
        (
            "operatingConditions",
            this->time().constant(),
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    cv_(transportProperties_.lookup("cv")),
    rho_(transportProperties_.lookup("rho")),
    K0_(transportProperties_.lookup("K")),
    k_(nullptr),
    sigma_(nullptr),
    DO2_(transportProperties_.lookup("DO2")),
    DV_(transportProperties_.lookup("DV")),
    epsilonP_(materialProperties_.lookup("epsilonP")),
    tau_(materialProperties_.lookup("tau")),
    MW_(materialProperties_.lookup("MW")),
    VW_(materialProperties_.lookup("VW")),
    sIm_(materialProperties_.lookup("sIm")),
    aLG_(materialProperties_.lookup("aLG")),
    p_(operatingConditions_.lookup("p")),
    RH_(operatingConditions_.lookup("RH")),
    c_
    (
        IOobject
        (
            "c",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this,
        dimensionedScalar("c0", dimensionSet(0, -3, 0, 0, 1, 0, 0), 1.0)
    ),
    DEffO2_(nullptr),
    DEffV_(nullptr),
    sRed_
    (
        IOobject
        (
            "sRed",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this,
        dimensionedScalar("sRed0", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.12)
    ),
    xVSat_
    (
        IOobject
        (
            "xVSat",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this,
        dimensionedScalar("xVSat0", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.1)
    ),
    K_(nullptr), 
    mu_
    (
        IOobject
        (
            "mu",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this,
        dimensionedScalar("mu", dimensionSet(0, 2, -1, 0, 0, 0, 0), 1.0)
    ), 
    dpCds_
    (
        IOobject
        (
            "dpCds",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this,
        dimensionedScalar("dpCds0", dimensionSet(1, -1, -2, 0, 0, 0, 0), 1.0)
    ),
    gamma_
    (
        IOobject
        (
            "gamma",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this,
        dimensionedScalar("gamma0", dimensionSet(0, 0, -1, 0, 0, 0, 0), 1.0)
    ), 
    sT_
    (
        IOobject
        (
            "sT",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this,
        dimensionedScalar("sT0", dimensionSet(1, -1, -3, 0, 0, 0, 0), 1.0)
    ),
    sV_
    (
        IOobject
        (
            "sV",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this,
        dimensionedScalar("sV0", dimensionSet(0, -3, -1, 0, 1, 0, 0), 1.0)
    ),
    ss_
    (
        IOobject
        (
            "ss",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this,
        dimensionedScalar("ss0", dimensionSet(0, -3, -1, 0, 1, 0, 0), 1.0)
    ),
    T_(nullptr),
    phiE_(nullptr),
    xO2_(nullptr),
    xV_(nullptr),
    s_(nullptr)
{
    k_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "k",
                this->time().timeName(),
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            *this,
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
                this->time().timeName(),
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            *this,
            dimensionedScalar(transportProperties_.lookup("sigma"))
        )
    );

    DEffO2_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "DEffO2",
                this->time().timeName(),
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            *this,
            DO2_
        )
    );

    DEffV_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "DEffV",
                this->time().timeName(),
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            *this,
            DV_
        )
    );

    K_.reset
    (
        new volScalarField
        (
            IOobject	
            (
                "K",
                this->time().timeName(),
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            *this,
            K0_
        )
    );

    T_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "T",
                this->time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        )
    );

    phiE_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "phiE",
                this->time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        )
    );

    xO2_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "xO2",
                this->time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        )
    );

    xV_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "xV",
                this->time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        )
    );

    s_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "s",
                this->time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        )
    );



    // thermal conductivity
    k_() = dimensionedScalar(transportProperties_.lookup("k"));
    // electric conducivity
    sigma_() = dimensionedScalar(transportProperties_.lookup("sigma"));
    
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::gasDiffusionLayer::~gasDiffusionLayer()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::gasDiffusionLayer::correct()
{
    // do nothing, add as required
}


void Foam::regionTypes::gasDiffusionLayer::setRDeltaT()
{
    // do nothing, add as required
}


void Foam::regionTypes::gasDiffusionLayer::setCoupledEqns()
{
    // update fields
    // gas species transport
    // ideal gas concentration
    c_ = p_/(RGas_*T_());
    // effective diffusion coefficient oxygen and vapor
    DEffO2_() = epsilonP_/sqr(tau_)*pow((1-s_()),3)*DO2_*pow((T_()/TRef_),1.5)*(pRef_/p_);
    DEffV_() = epsilonP_/sqr(tau_)*pow((1-s_()),3)*DV_*pow((T_()/TRef_),1.5)*(pRef_/p_);
    
    // liquid water transport
    // reduced liquid water saturation
    sRed_ = (s_()-sIm_)/(1-sIm_);
    // saturation vapor fraction
    xVSat_ = (exp(23.1963-(TRefP1_/(T_()-TRefP2_)))*pDim_)/p_;
    // dynamic viscosity water
    mu_ = exp(-3.63148+(TRefMu1_/(T_()+TRefMu2_)))*muDim_;
    // derivate of capillary pressure with respect to liquid water saturation
    dpCds_ = (4.8422e-3*exp(-44.02*(s_()-0.496))+2255.0649*exp(8.103*(s_()-0.496)))*pDim_;
    // reduced liquid water permeability
    K_() = (1e-6+pow(sRed_,3))*K0_;
    // evaporation/condensation rate
    {if(xV_() < xVSat_) //evaporation
	{
	   gamma_ = 5e-4*sqrt(RGas_*T_()/(2*pi_*MW_))*aLG_*sRed_;
	}
	else //condensation
	{
	   gamma_ = 6e-3*sqrt(RGas_*T_()/(2*pi_*MW_))*aLG_*sRed_;
	}}

    // source terms
    // heat Source - joule heating electrons and condensation/evaporation heat
    //sT_ = sigma_()*magSqr(fvc::grad(phiE_()))+gamma_*c_*(xV_()-xVSat_)*HEC_;
    // mass source vapor
    //sV_ = -gamma_*c_*(xV_()-xVSat_);
    // mass source liquid water
    //ss_ = gamma_*c_*(xV_()-xVSat_);
    // fourier heat conduction
    fvScalarMatrix TEqn =
    (
        rho_*cv_*fvm::ddt(T_())
     ==
        fvm::laplacian(k_(), T_(), "laplacian(k,T)")
	//+sT_
    );

    // set Eqns
    // ohm's law for electrons
    fvScalarMatrix phiEEqn =
    (
        -fvm::laplacian(sigma_(), phiE_(), "laplacian(sigma,phiE)")
    );

    // fick diffusion for oxygen
    fvScalarMatrix xO2Eqn =
    (
	c_*fvm::ddt(xO2_())
     ==
        fvm::laplacian(c_*DEffO2_(), xO2_(), "laplacian(D,x)")
    );

    // fick diffusion for vapor
    fvScalarMatrix xVEqn =
    (
	c_*fvm::ddt(xV_())
     ==
        fvm::laplacian(c_*DEffV_(), xV_(), "laplacian(D,x)")
       //+sV_
    );

    // liquid water transport (derived from Darcy's Law)
    fvScalarMatrix sEqn =
    (
	1/MW_*fvm::ddt(s_())
     ==
        fvm::laplacian(K_()/(mu_*VW_)*dpCds_, s_(), "laplacian(K,s)")
       //+ss_
    );

    fvScalarMatrices.set
    (
        T_().name() + this->name() + "Eqn",
	new fvScalarMatrix(TEqn)
    );

    fvScalarMatrices.set
    (
	phiE_().name() + this->name() + "Eqn",
	new fvScalarMatrix(phiEEqn)
    );

    fvScalarMatrices.set
    (
	xO2_().name() + this->name() + "Eqn",
	new fvScalarMatrix(xO2Eqn)
    );

    fvScalarMatrices.set
    (
	xV_().name() + this->name() + "Eqn",
	new fvScalarMatrix(xVEqn)
    );

    fvScalarMatrices.set
    (
	s_().name() + this->name() + "Eqn",
	new fvScalarMatrix(sEqn)
    );
}

void Foam::regionTypes::gasDiffusionLayer::updateFields()
{
    // do nothing, add as required
}

void Foam::regionTypes::gasDiffusionLayer::solveRegion()
{
    // do nothing, add as required
}


// ************************************************************************* //
