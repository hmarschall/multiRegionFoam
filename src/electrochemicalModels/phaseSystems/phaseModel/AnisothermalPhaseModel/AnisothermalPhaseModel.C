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

#include "AnisothermalPhaseModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::AnisothermalPhaseModel<BasePhaseModel>::AnisothermalPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, index),
    kappa_
    (
            IOobject
            (
                IOobject::groupName("kappa", this->name()),
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::READ_IF_PRESENT,        // I need kappa also updated at the boundary
                IOobject::AUTO_WRITE        // Dunno how to set the zeroGrad BC otherwise for now
             ),
            this->mesh(),
            dimensionedScalar("kappa", dimensionSet(1, 1, -3,-1,0,0,0), 1)
    )
    // kappa_(nullptr)
    // T_(nullptr)
{
    // kappa_.reset
    // (
    //     new volScalarField
    //     (
    //         IOobject
    //         (
    //             IOobject::groupName("kappa", this->name()),
    //             this->mesh().time().timeName(),
    //             this->mesh(),
    //             IOobject::READ_IF_PRESENT,        // I need kappa also updated at the boundary
    //             IOobject::AUTO_WRITE        // Dunno how to set the zeroGrad BC otherwise for now
    //          ),
    //         this->mesh(),
    //         dimensionedScalar("kappa", dimensionSet(1, 1, -3,-1,0,0,0), 1)
    //     )
    // );

    // T_.reset
    // (
    //     new volScalarField
    //     (
    //         IOobject
    //         (
    //             IOobject::groupName("T", this->name()),
    //             this->mesh().time().timeName(),
    //             this->mesh(),
    //             IOobject::MUST_READ,        // I need kappa also updated at the boundary
    //             IOobject::AUTO_WRITE        // Dunno how to set the zeroGrad BC otherwise for now
    //         ),
    //         this->mesh()
    //     )
    // );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::AnisothermalPhaseModel<BasePhaseModel>::~AnisothermalPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::AnisothermalPhaseModel<BasePhaseModel>::correctThermo()
{
    BasePhaseModel::correctThermo();

    // kappa_.correctBoundaryCondiitons(); // So that there are 

    volScalarField& T = this->thermo_->T();
    T.correctBoundaryConditions();

    this->thermo_->correct();
}


template<class BasePhaseModel>
bool Foam::AnisothermalPhaseModel<BasePhaseModel>::isothermal() const
{
    return false;
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::AnisothermalPhaseModel<BasePhaseModel>::heEqn()
{
    const volScalarField& alpha = *this;

    const volVectorField& U(this->U());
    const surfaceScalarField alphaPhi(this->alphaPhi());
    const surfaceScalarField alphaRhoPhi(this->alphaRhoPhi());

    const volScalarField contErr(this->continuityError());
    // const volScalarField K(this->K());

    // volScalarField& he = this->thermo_->he();
    volScalarField& he = this->thermo_->he();    // update thermo model

    //- Porosity
    //- Non-dimension, created via alpha
    volScalarField porosity
    (
        IOobject
        (
            "porosity",
            this->fluid().mesh().time().timeName(),
            this->fluid().mesh()
        ),
        this->fluid().mesh(),
        dimensionedScalar("por", dimless, 1.0)
    );

    // if (!this->fluid().porosityModels().empty())
    // {
    //     forAll(this->fluid().porosityModels(), iz)
    //     {
    //         label znId = this->fluid().mesh().
    //             cellZones().findZoneID(this->fluid().porosityModels()[iz].zoneName());

    //         scalar por = this->fluid().porosityModels()[iz].porosity();

    //         labelList znCells(this->fluid().mesh().cellZones()[znId]);

    //         forAll(znCells, cellI)
    //         {
    //             label cellId = znCells[cellI];

    //             porosity[cellId] *= por;
    //         }
    //     }
    // }


    tmp<fvScalarMatrix> tEEqn
    (
        fvm::ddt(alpha*this->thermo_->rho(), he)
      + fvm::div(alphaRhoPhi, he)
      + fvm::SuSp(-contErr, he)
    //   - contErr*K
      - fvm::laplacian
        (
            // fvc::interpolate(alpha*porosity)
            fvc::interpolate(alpha)
           *fvc::interpolate(this->alphaEff()), // missing, it is alpha turbulent part
            he
        )
    //    +  this->thermophysicalTransport_->divq(he)
        // alpha*this->Qdot()      // TODO: This would be reasonable if I do not have a main region
    );

    Info << "Debug line anisothermalPhaseModel " << endl;

    return tEEqn;
}

template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::AnisothermalPhaseModel<BasePhaseModel>::TEqn()
{
    const volScalarField& alpha = *this;

    volScalarField& T = this->thermo_->T();
    // volScalarField& he = this->thermo_->he();

    //- Porosity
    //- Non-dimension, created via alpha
    // TODO: For now take the laminar heat conductivity 
    // == thermal diffusivity for temperature
    // volScalarField& kappaEff = this->thermo_->kappa()();
    // kappa_() = this->thermo_->kappa()();
    tmp<volScalarField> tkappa = this->thermo_->kappa();
    kappa_ = tkappa();

    tmp<volScalarField> tcpv = this->thermo_->Cpv();
    volScalarField& cpv = tcpv();

    const volScalarField& rho = this->thermo().rho();

    const volVectorField& U = this->URef();
    const tmp<surfaceScalarField> alphaPhi(this->alphaPhi());
    const tmp<surfaceScalarField> alphaRhoPhi(this->alphaRhoPhi());
    // const surfaceScalarField alphaRhoCpPhi(alphaRhoPhi*linearInterpolate(cpv));
    const surfaceScalarField faceCPV(linearInterpolate(cpv));
    const surfaceScalarField alphaRhoCpPhi(alphaRhoPhi()*faceCPV);

    Info << "Print average(alphaRhoPhi()) " << average(alphaRhoPhi()) << endl;
    Info << "Print average(alphaRhoCpPhi) " << average(alphaRhoCpPhi) << endl;

        //- Porosity
    //- Non-dimension, created via alpha
    volScalarField porosity
    (
        IOobject
        (
            "porosity",
            this->fluid().mesh().time().timeName(),
            this->fluid().mesh()
        ),
        this->fluid().mesh(),
        dimensionedScalar("por", dimless, 1.0)
    );


    // if (!this->fluid().porosityModels().empty())
    // {
    //     forAll(this->fluid().porosityModels(), iz)
    //     {
    //         label znId = this->mesh().
    //             cellZones().findZoneID(this->fluid().porosityModels()[iz].zoneName());

    //         scalar por, kZn;
    //         this->fluid().porosityModels()[iz].dict().lookup("porosity") >> por;
    //         this->fluid().porosityModels()[iz].dict().lookup("k") >> kZn;

    //         labelList znCells(this->mesh().cellZones()[znId]);

    //         forAll(znCells, cellI)
    //         {
    //             label cellId = znCells[cellI];

    //             kappa_[cellId] *= por;
    //             kappa_[cellId] += kZn*(1.0 - por)/cpv[cellId];
    //         }
    //     }
    // }

    Info << "Print alpha " << average(alpha) << endl;

    tmp<fvScalarMatrix> tEEqn
    (
        fvm::ddt(alpha*rho*cpv, T)
      + fvm::div(alphaRhoCpPhi, T, "div(rhoCpPhi,T)")
      - fvm::laplacian
        (
            fvc::interpolate(alpha*porosity)
            *fvc::interpolate(kappa_),
            T
        )
     ==
        alpha*this->Qdot()
    );

    Info << "Heat flux due to the electrochemical reactions " << average(this->Qdot()()) << endl; 

    return tEEqn;
}

template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::AnisothermalPhaseModel<BasePhaseModel>::heQdot()
{
    const volScalarField& alpha = *this;

    const volVectorField U(this->U());
    const surfaceScalarField alphaPhi(this->alphaPhi());
    const surfaceScalarField alphaRhoPhi(this->alphaRhoPhi());

    const volScalarField contErr(this->continuityError());
    const volScalarField K(this->K());

    // For foam extend
    const volScalarField alphaRho = alpha*this->rho();

    volScalarField& he = this->thermo_->he();

    tmp<volScalarField> qdot
    (
    //   - fvc::ddt(alpha, this->rho(), K) 
      - fvc::ddt(alphaRho, K) 
      - fvc::div(alphaRhoPhi, K)
      + contErr*K
//
//      + alpha*this->Qdot()
    );

    // // Add the appropriate pressure-work term
    // if (he.name() == this->thermo_->phasePropertyName("e"))
    // {
    //     qdot.ref() -= filterPressureWork
    //     (
    //         fvc::div(fvc::absolute(alphaPhi, alpha, U), this->thermo().p())
    //       + this->thermo().p()*fvc::ddt(alpha)
    //     );
    // }
    // else if (this->thermo_->dpdt())
    // {
    //     qdot.ref() += filterPressureWork(alpha*this->fluid().dpdt());
    // }

    return qdot;
}


// ************************************************************************* //
