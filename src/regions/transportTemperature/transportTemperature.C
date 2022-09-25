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
#include "transportTemperature.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(transportTemperature, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        transportTemperature,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::transportTemperature::transportTemperature
(
    const Time& runTime,
    const word& regionName
)
:
    regionType(runTime, regionName),

    regionName_(regionName),

    // U_
    // (
    //     IOobject
    //     (
    //         "U",
    //         mesh().time().timeName(),
    //         mesh(),
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh()
    // ),
//    phi_
//    (
//		IOobject
//		(
//			"phi",
//            mesh().time().timeName(),
//            mesh(),
//			IOobject::READ_IF_PRESENT,
//			IOobject::AUTO_WRITE
//		),
//		linearInterpolate(U_) & mesh().Sf()    
//    ),

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
    cp_(transportProperties_.lookup("cp")),
    rho_(transportProperties_.lookup("rho")),

//    alpha_
//    (
//        IOobject
//        (
//            "alpha",
//            this->time().timeName(),
//            *this,
//            IOobject::READ_IF_PRESENT,
////            IOobject::MUST_READ,
//            IOobject::NO_WRITE
//        ),
//        *this,
//        k_/(rho_*cp_)
//    ),
//    T_
//    (
//        IOobject
//        (
//            "T",
//            this->time().timeName(),
//            *this,
//            IOobject::MUST_READ,
//            IOobject::AUTO_WRITE
//        ),
//        *this
//    )
    U_(nullptr),
    k_(nullptr),
    alpha_(nullptr),
    phi_(nullptr),
    T_(nullptr)
{
    U_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "U",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh()
            )
        );

    phi_.reset
        (
            new surfaceScalarField
            (
		        IOobject
		        (
			        "phi",
                    mesh().time().timeName(),
                    mesh(),
			        IOobject::NO_READ,
			        IOobject::AUTO_WRITE
		        ),
		        linearInterpolate(U_()) & mesh().Sf()
            )
        );

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

    alpha_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "alpha",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("alpha", dimensionSet(0,2,-1,0,0,0,0), 0)
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

    alpha_() = k_()/(rho_*cp_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::transportTemperature::~transportTemperature()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::transportTemperature::correct()
{
    // do nothing, add as required
}


void Foam::regionTypes::transportTemperature::setRDeltaT()
{
    // do nothing, add as required
}


void Foam::regionTypes::transportTemperature::setCoupledEqns()
{
    fvScalarMatrix TEqn =
    (
        fvm::ddt(T_())
      + fvm::div(phi_(), T_())
     ==
        fvm::laplacian(alpha_(), T_())
    );

    fvScalarMatrices.set
    (
        T_().name() + mesh().name() + "Eqn",
        new fvScalarMatrix(TEqn)
    );
}

void Foam::regionTypes::transportTemperature::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportTemperature::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportTemperature::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportTemperature::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportTemperature::pressureCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //
