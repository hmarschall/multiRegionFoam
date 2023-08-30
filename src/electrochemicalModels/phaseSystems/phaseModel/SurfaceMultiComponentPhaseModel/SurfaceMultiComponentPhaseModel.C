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

#include "SurfaceMultiComponentPhaseModel.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvMatrix.H"

#include "zeroGradientFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::SurfaceMultiComponentPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, index),
    rhoField_
    (
        IOobject
        (
            "rhoField." + phaseName,
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->thermo().rho()
        // this->mesh(),
        // dimensionedScalar("rhoField", dimDensity, 0),
        // zeroGradientFvPatchScalarField::typeName
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        // readScalar(this->mesh().solverDict("Yi"))
        // readScalar(this->mesh().solutionDict("Yi").lookup("residualAlpha"))
        1e-6    // foam-extend does not now solverDict
    ),
    inertIndex_(-1)
    // thermophysicalTransport_
    // (
    //     ThermophysicalTransportModel<compressibleMomentumTransportModel,this->thermo_->thermoName()>::New
    //     (
    //         this->momentumTransport_,
    //         this->thermo_
    //     )
    // )
{

    // Create the multiSpecieModel
    thermoPhysicalTransportModel_ = multiSpeciesTransportModel::New(*this); 

    const word inertSpecie
    (
        this->thermo_->lookupOrDefault("inertSpecie", word::null)
    );

    if (inertSpecie != word::null)
    {
        inertIndex_ = this->thermo_->composition().species()[inertSpecie];
    }

    PtrList<volScalarField>& Y = this->thermo_->composition().Y();

    //- Mixture mole fraction
    const volScalarField W(this->thermo_->W());

    //- Mole fraction, X
    X_.resize(Y.size());
    iDmdt_.resize(Y.size());

    forAll(Y, i)
    {
        const dimensionedScalar Wi
        (
            "W",
            dimMass/dimMoles,
            this->thermo_->composition().Wi(i)
        );

        X_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("X", Y[i].member()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                W*Y[i]/Wi
            )
        );

        iDmdt_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("iDmdt", Y[i].member()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("iDmdt", dimDensity/dimTime, 0.0)
            )
        );

        if (i != inertIndex_)
        {
            const label j = YActive_.size();
            YActive_.resize(j + 1);
            YActive_.set(j, &Y[i]);

            XActive_.resize(j + 1);
            XActive_.set(j, &X_[i]);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::~SurfaceMultiComponentPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::correctThermo()
{
  BasePhaseModel::correctThermo();

    rhoField_ = this->thermo_->rho();
    rhoField_.correctBoundaryConditions();

    volScalarField Yt
    (
        IOobject
        (
            IOobject::groupName("Yt", this->name()),
            this->mesh().time().timeName(),
            this->mesh()
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 0)
    );

    PtrList<volScalarField>& Yi = YRef();

    forAll(Yi, i)
    {
        if (i != inertIndex_)
        {
            Yt += Yi[i];
        }
    }

    if (inertIndex_ != -1)
    {
        Yi[inertIndex_] = scalar(1) - Yt;
        Yi[inertIndex_].max(0);

        Yt += Yi[inertIndex_];
    }

    //- Normalize
    forAll(Yi, i)
    {
        Yi[i] /= Yt;

        Info<< "Y: " << Yi[i].name() << ":"
            << gMin(Yi[i]) << " -> " << gMax(Yi[i])
            << endl;
    }

    //- Mixture mole fraction
    const volScalarField W(this->thermo_->W());

    //- Update mole fraction X_
    forAll(Yi, i)
    {
        const dimensionedScalar Wi
        (
            "W",
            dimMass/dimMoles,
            this->thermo_->composition().Wi(i)
        );

        X_[i] = W*Yi[i]/Wi;
    }

    thermoPhysicalTransportModel_->correct();
}


template<class BasePhaseModel>
void Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::correctBC()
{
    // Calculate species sources and sinks and set
    // boundary conditions for mass fractions and velocity
    // at fluid/electrode interfaces

    // For the air side - is this needed!? -> I think not really

    Info << "Correct the fixedGradient boundary conditions"
         <<  " for the species mass flux (multicomponentPhaseModel) " << endl;

    const dictionary& dict = this->fluid(); 
    const word fluidCatalystName = word(dict.lookup("fluidCatalystName")); 
    label fluidCatalystID = this->mesh().boundaryMesh().findPatchID(fluidCatalystName);

    PtrList<volScalarField>& Y = this->thermo_->composition().Y();
    const volScalarField& rho = this->thermo().rho();

    // Calculate the BC for the mass fraction of each species except the inert specie 
    forAll(Y, i)
    {
        label speciesIndex = this->thermo_->composition().species()[Y[i].member()];

        if(speciesIndex != inertIndex_)
        {
            // The dissolved water flux is already integrated into iDmdt(H2O)
            const volScalarField& idmdt = iDmdt(Y[i].member());
            const fvPatchField<scalar>& massFluxI = idmdt.boundaryField()[fluidCatalystID];
            tmp<volScalarField> tDiffRho = thermoPhysicalTransportModel_->DEff(Y[i]);
            volScalarField& diffRho = tDiffRho();

            const scalarField diffRhoPatch =
            ( 
                diffRho.boundaryField()[fluidCatalystID]
            );

            volScalarField& Yi = Y[i];
            fixedGradientFvPatchScalarField& YiBC =
            refCast<fixedGradientFvPatchScalarField>
            (
                Yi.boundaryField()[fluidCatalystID]
            );

            YiBC.gradient() = massFluxI *(1.0 - YiBC) ;

            // Add changes due to other species
            forAll(Y, j)
            {
               if((j != i)) //  && isFlux[j])
               {
                    const volScalarField& idmdtj = iDmdt(Y[j].member());
                    const fvPatchField<scalar>& massFluxJ = idmdtj.boundaryField()[fluidCatalystID];

                    YiBC.gradient() -= massFluxJ*YiBC;
               }
            }

            YiBC.gradient() /= (diffRhoPatch) ; 
        }
    }

    // // Reference to the phase velocity
    // volVectorField& U = this->URef();

    // const scalarField& rhoPatch = rho.boundaryField()[fluidCatalystID];

    // // Summation of mass fluxes of species within this phase
    // volScalarField mfluxSum = iDmdt_[0] * 0;
    // fvPatchField<scalar>& mfluxSumPatch = mfluxSum.boundaryField()[fluidCatalystID];

    // forAll(Y, i)
    // {
    //     mfluxSumPatch += iDmdt(Y[i].member()).boundaryField()[fluidCatalystID];
    // }

    // U.boundaryField()[fluidCatalystID] ==
    // (
    //     -(mfluxSumPatch/rhoPatch)
    //     *(this->mesh().Sf().boundaryField()[fluidCatalystID])
    //     /(this->mesh().magSf().boundaryField()[fluidCatalystID])
    // ); 
}



template<class BasePhaseModel>
bool Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::pure() const
{
    return false;
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::YiEqn(volScalarField& Yi)
{

    const volScalarField& alpha = *this;
    const surfaceScalarField& alphaRhoPhi(this->alphaRhoPhi());
    const volScalarField& rho = this->thermo().rho();

    const volScalarField alphaRho = alpha*rho;

    label speciesIndex = this->thermo_->composition().species()[Yi.member()];

    Info << "Print speciesIndex " << speciesIndex << endl;
    Info << "Print name of the specie " << Yi.member() << endl;
    Info << "Print average(iDmdt_[speciesIndex]) " << average(iDmdt_[speciesIndex]) << endl;

    return
    (
        // fvm::ddt(alpha, rho, Yi)
        fvm::ddt(alphaRho, Yi)
      + fvm::div(alphaRhoPhi, Yi, "div(" + alphaRhoPhi.name() + ",Yi)")
      + thermoPhysicalTransportModel_->divj(Yi)
    //   ==

    //     fvc::ddt(residualAlpha_*rho, Yi)
    //   - fvm::ddt(residualAlpha_*rho, Yi)
    );
}


template<class BasePhaseModel>
const Foam::PtrList<Foam::volScalarField>&
Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::Y() const
{
    return this->thermo_->composition().Y();
}


template<class BasePhaseModel>
const Foam::PtrList<Foam::volScalarField>&
Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::iDmdt() const
{
    return iDmdt_;
}


template<class BasePhaseModel>
const Foam::PtrList<Foam::volScalarField>&
Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::X() const
{
    return X_;
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::Y(const word& name) const
{
    return this->thermo_->composition().Y(name);
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::iDmdt(const word& name) const
{
    return iDmdt_[this->thermo_->composition().species()[name]];
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::X(const word& name) const
{
    return X_[this->thermo_->composition().species()[name]];
}


template<class BasePhaseModel>
Foam::PtrList<Foam::volScalarField>&
Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::YRef()
{
    return this->thermo_->composition().Y();
}


template<class BasePhaseModel>
Foam::PtrList<Foam::volScalarField>&
Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::iDmdtRef()
{
    return iDmdt_;
}


template<class BasePhaseModel>
Foam::PtrList<Foam::volScalarField>&
Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::XRef()
{
    return X_;
}


template<class BasePhaseModel>
const Foam::UPtrList<Foam::volScalarField>&
Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::YActive() const
{
    return YActive_;
}


template<class BasePhaseModel>
const Foam::UPtrList<Foam::volScalarField>&
Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::XActive() const
{
    return XActive_;
}


template<class BasePhaseModel>
Foam::UPtrList<Foam::volScalarField>&
Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::YActiveRef()
{
    return YActive_;
}


template<class BasePhaseModel>
Foam::UPtrList<Foam::volScalarField>&
Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::XActiveRef()
{
    return XActive_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::SurfaceMultiComponentPhaseModel<BasePhaseModel>::dmdt() const
{
    tmp<volScalarField> dmdt
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("dmdt0", this->name()),
                this->mesh().time().timeName(),
                this->mesh()
            ),
            this->mesh(),
            dimensionedScalar("zero", dimDensity/dimTime, 0)
        )
    );

    scalarField& dmdt0 = dmdt();

    forAll(Y(), i)
    {
        volScalarField& Yi = const_cast<volScalarField&>(Y()[i]);

        dmdt0 += iDmdt_[i];
    }

    return dmdt;
}



// ************************************************************************* //
