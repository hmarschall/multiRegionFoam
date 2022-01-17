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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::ionomer::ionomer
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
	    this->time().timeName(),
	    *this,
	    IOobject::READ_IF_PRESENT,
	    IOobject::NO_WRITE
	),
	*this,
	dimensionedScalar("f0", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.06)
    ),
    fCond_
    (
	IOobject
	(
	    "fCond",
	    this->time().timeName(),
	    *this,
	    IOobject::READ_IF_PRESENT,
	    IOobject::NO_WRITE
	),
	*this,
	dimensionedScalar("f0", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.06)
    ),
    xi_
    (
	IOobject
	(
	    "xi",
	    this->time().timeName(),
	    *this,
	    IOobject::READ_IF_PRESENT,
	    IOobject::NO_WRITE
	),
	*this,
	dimensionedScalar("xi0", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0)
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
	dimensionedScalar("sT0", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
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
                this->time().timeName(),
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            *this,
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
                this->time().timeName(),
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            *this,
            kappa0_
        )
    );

    DLambda_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "DLambda",
                this->time().timeName(),
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            *this,
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
                this->time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        )
    );

    phiP_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "phiP",
                this->time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        )
    );

    lambda_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "lambda",
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

    // update variables

    // ionomer properties
    // water diffusion through ionomer
    DLambda_() = (3.842*pow(lambda_(),3)-32.03*sqr(lambda_())+67.75*lambda_())/(pow(lambda_(),3)-2.115*sqr(lambda_())-33.013*lambda_()+103.37)*DLambda0_*exp(ELambda_/RGas_*((1/TRef_)-(1/T_())));
    // electro-osmotic drag coefficient
    xi_ = 2.5*lambda_()/22;
    // volume fraction water and protonic conductivity
    f_ = lambda_()*VW_/(lambda_()*VW_+VM_);
    fCond_ = fComp_;
    {if(f_ > fCond_)
    {
	kappa_() = kappa0_*pow((f_-0.06),1.5)*exp(EKappa_/RGas_*((1/TRef_)-(1/T_())));
    }
    else
    {
	kappa_() = kappa0_*pow(0,1.5)*exp(EKappa_/RGas_*((1/TRef_)-(1/T_())));
    }}

    // source terms
    // heat source - joule heating protons
    //sT_ = kappa_()*magSqr(fvc::grad(phiP_()));
    
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
        fvm::laplacian(DLambda_()/VM_, lambda_(), "laplacian(DLambda,lambda)")+fvm::laplacian(xi_*kappa_()/FConst_, phiP_(), "laplacian(kappa,phiP)")
    );
    
    fvScalarMatrices.set
    (
        T_().name() + this->name() + "Eqn",
	new fvScalarMatrix(TEqn)
    );

    fvScalarMatrices.set
    (
	phiP_().name() + this->name() + "Eqn",
	new fvScalarMatrix(phiPEqn)
    );

    fvScalarMatrices.set
    (
	lambda_().name() + this->name() + "Eqn",
	new fvScalarMatrix(lambdaEqn)
    );
}

void Foam::regionTypes::ionomer::updateFields()
{
    // do nothing, add as required
}

void Foam::regionTypes::ionomer::solveRegion()
{
    // do nothing, add as required
}


// ************************************************************************* //
