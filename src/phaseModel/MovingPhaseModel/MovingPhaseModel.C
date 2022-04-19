/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "MovingPhaseModel.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::phi(const volVectorField& U) const
{
    word phiName("phi");

    IOobject phiHeader
    (
        phiName,
        U.mesh().time().timeName(),
        U.mesh(),
        IOobject::NO_READ
    );

    if (phiHeader.headerOk())
    {
        Info<< "Reading face flux field " << phiName << endl;

        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    U.mesh().time().timeName(),
                    U.mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                U.mesh()
            )
        );
    }
    else
    {
        Info<< "Calculating face flux field " << phiName << endl;

        wordList phiTypes
        (
            U.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        forAll(U.boundaryField(), i)
        {
            if
            (
                isA<fixedValueFvPatchVectorField>(U.boundaryField()[i])
             || isA<slipFvPatchVectorField>(U.boundaryField()[i])
             || isA<partialSlipFvPatchVectorField>(U.boundaryField()[i])
            )
            {
                phiTypes[i] = fixedValueFvPatchScalarField::typeName;
            }
        }

        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    U.mesh().time().timeName(),
                    U.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(U) & U.mesh().Sf(),
                phiTypes
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MovingPhaseModel<BasePhaseModel>::MovingPhaseModel
(
        const dictionary& dict,
        const fvMesh& mesh
)
:
    BasePhaseModel(dict, mesh),
    U_
    (
        IOobject
        (
            "U",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    phi_(phi(U_)),
    alphaPhi_
    (
        IOobject
        (
            "alphaPhi",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0)
    ),
    alphaRhoPhi_
    (
        IOobject
        (
            "alphaRhoPhi",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(1, 0, -1, 0, 0), 0)
    ),
    divU_(nullptr),
    continuityErrorFlow_
    (
        IOobject
        (
            "continuityErrorFlow",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimDensity/dimTime, 0)
    ),
    continuityErrorSources_
    (
        IOobject
        (
            "continuityErrorSources",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimDensity/dimTime, 0)
    ),
    K_(nullptr)
{
    phi_.writeOpt() = IOobject::AUTO_WRITE;

    correctKinematics();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MovingPhaseModel<BasePhaseModel>::~MovingPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::correct()
{
    BasePhaseModel::correct();

    volScalarField& rho = this->thermoRef().rho();

    continuityErrorFlow_ =  fvc::div(alphaRhoPhi_);
}


template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::correctKinematics()
{
    BasePhaseModel::correctKinematics();

    if (K_.valid())
    {
        K_.clear();
        K();
    }
}


template<class BasePhaseModel>
bool Foam::MovingPhaseModel<BasePhaseModel>::stationary() const
{
    return false;
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::MovingPhaseModel<BasePhaseModel>::UEqn()
{
    const volScalarField& alpha = *this;
    const volScalarField& rho = this->thermo().rho();

    // tmp<fvVectorMatrix> UEqn
    // (
    //     fvm::div(alphaRhoPhi_, U_)
    //   + fvm::SuSp(- this->continuityError(), U_)
    // );

    // // Introduce darcy term here (different contribution for each porous cell zone)
    // this->porosityModels().addResistance(UEqn);

    // return
    // (
    //    UEqn
    // );

    return
    (
        fvm::div(alphaRhoPhi_, U_)
      + fvm::SuSp(- this->continuityError(), U_)
    );
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::MovingPhaseModel<BasePhaseModel>::UfEqn()
{
    // As the "normal" U-eqn but without the ddt terms

    const volScalarField& alpha = *this;
    const volScalarField& rho = this->thermo().rho();

    // tmp<fvVectorMatrix> UEqn
    // (
    //     fvm::div(alphaRhoPhi_, U_)
    //   - fvm::Sp(fvc::div(alphaRhoPhi_), U_)
    //   + fvm::SuSp(- this->continuityErrorSources(), U_)

    // );

    // // Introduce darcy term here
    // this->porosityModels().addResistance(UEqn);

    // return
    // (
    //    UEqn
    // );

    return
    (
        fvm::div(alphaRhoPhi_, U_)
      - fvm::Sp(fvc::div(alphaRhoPhi_), U_)
      + fvm::SuSp(- this->continuityErrorSources(), U_)
    );
}


template<class BasePhaseModel>
Foam::tmp<Foam::volVectorField>
Foam::MovingPhaseModel<BasePhaseModel>::U() const
{
    return U_;
}


template<class BasePhaseModel>
Foam::volVectorField&
Foam::MovingPhaseModel<BasePhaseModel>::URef()
{
    return U_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::phi() const
{
    return phi_;
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::MovingPhaseModel<BasePhaseModel>::phiRef()
{
    return phi_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::alphaPhi() const
{
    return alphaPhi_;
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::MovingPhaseModel<BasePhaseModel>::alphaPhiRef()
{
    return alphaPhi_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::alphaRhoPhi() const
{
    return alphaRhoPhi_;
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::MovingPhaseModel<BasePhaseModel>::alphaRhoPhiRef()
{
    return alphaRhoPhi_;
}

template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::continuityError() const
{
    return continuityErrorFlow_ + continuityErrorSources_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::continuityErrorFlow() const
{
    return continuityErrorFlow_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::continuityErrorSources() const
{
    return continuityErrorSources_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::K() const
{
    if (!K_.valid())
    {
        volScalarField& K0 = K_();
        // K_ =
        //     new volScalarField
        //     (
        //         // IOobject
        //         // (
        //         //     "K",
        //         //     this->mesh().time().timeName(),
        //         //     this->mesh(),
        //         //     IOobject::NO_READ,
        //         //      IOobject::NO_WRITE
        //         // ),
        //         "K",
        //         0.5*magSqr(this->U())
        //     );
        K0 = 0.5*magSqr(this->U());
    }

    return tmp<volScalarField>(K_());
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::divU() const
{
    return divU_.valid() ? tmp<volScalarField>(divU_()) : tmp<volScalarField>();
}


template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::divU(tmp<volScalarField> divU)
{
    divU_ = divU;
}

// ************************************************************************* //
