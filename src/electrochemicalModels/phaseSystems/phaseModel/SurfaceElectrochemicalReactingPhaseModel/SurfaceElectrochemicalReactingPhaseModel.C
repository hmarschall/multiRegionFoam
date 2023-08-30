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

#include "SurfaceElectrochemicalReactingPhaseModel.H"
#include "phaseSystem.H"
#include "surfaceDissolvedModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::SurfaceElectrochemicalReactingPhaseModel<BasePhaseModel>::SurfaceElectrochemicalReactingPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, index),
    eta_
    (
        IOobject
        (
            "eta",
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar
        (
            "eta",
            dimensionSet(1, 2, -3, 0, 0, -1, 0),
            0.0
        )
    ),
    jump_
    (
        IOobject
        (
            // IOobject::groupName("jump", phase.mesh().name()),
            "jump",
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar
        (
            "jump",
            dimensionSet(1, 2, -3, 0, 0, -1, 0),
            0.0
        )
    ),
    regions_(fluid.subDict("regions")),
    saturation_(saturationModel::New(fluid.subDict("saturation"), fluid.mesh())),
    phiNames_(fluid.lookup("phiNames")),
    iNames_(fluid.lookup("iNames")),
    relax_(fluid.lookupOrDefault<scalar>("relax", 1.0)),
    alpha_(fluid.lookupOrDefault<scalar>("alpha", 0.5)),
    gamma_(fluid.lookupOrDefault<scalar>("gamma", 1.0)),
    i0Value_("i0", dimCurrent/dimArea, readScalar(fluid.lookup("i0"))), 
    i_
    (
        IOobject
        (
            IOobject::groupName("i", fluid.mesh().name()),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar
        (
            "i",
            dimCurrent/dimArea,
            0.0
        )
    ),
    i0_
    (
        IOobject
        (
            "i0",
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar
        (
            "i0",
            dimCurrent/dimArea,
            0.0
        )
    ),
    rxnList_(fluid.lookup("RxnList")),
    // Nernst variables
    nernst_
    (
        IOobject
        (
            "nernst",
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar
        (
            "nernst",
            dimensionSet(1, 2, -3, 0, 0, -1, 0),
            0.0
        )
    ),
    deltaH_
    (
        IOobject
        (
            "deltaH",
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar
        (
            "deltaH",
            dimEnergy/dimMoles,
            0.0
        )
    ),
    deltaS_
    (
        IOobject
        (
            "deltaS",
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar
        (
            "deltaS",
            dimEnergy/dimMoles/dimTemperature,
            0.0
        )
    ),
    pRef_("pRef", dimPressure, readScalar(fluid.lookup("pRef"))),
    species_(fluid.subDict("species")),
    A_(0),
    B_(0),
    Qdot_
    (
        IOobject
        (
            IOobject::groupName("Qdot", this->name()),
            eta_.mesh().time().timeName(),
            eta_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        eta_.mesh(),
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0)
    )
{
    // Info << "Breakpoint SurfaceElectrochemicalReactingPhaseModel.C " << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::SurfaceElectrochemicalReactingPhaseModel<BasePhaseModel>::~SurfaceElectrochemicalReactingPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::SurfaceElectrochemicalReactingPhaseModel<BasePhaseModel>::correct()
{
    //- Activation overpotential model update
    correctOverpotential();

    const PtrList<volScalarField>& Y = this->thermo_->composition().Y();
    
    word water(phaseModel::water);

    forAll(Y,i)
    {
        if (Y[i].member() != water)
        {
            const volScalarField& Yi = Y[i];
            // Correct the dmdt[i] of the gases
            correctGasIdmdt(Yi);
        }

    }

    //- Get sub regions
    //- Refer to regionType
    //- Including: fluid, electron (BPP + GDL + CL), ion (CLs + membrane)
    const regionType& fluidPhase = this->region
    (
        word(this->regions_.subDict("fluid").lookup("name"))
    );

    const regionType& ionPhase = this->region
    (
        word(this->regions_.subDict("ion").lookup("name"))
    );

    // // Get the phase System
    //  const phaseSystem& phaseSys = this->fluid().mesh().template
    //  lookupObject<phaseSystem>(phaseSystem::propertiesName);

    // //- Get the present phase Model from fluid phase
    // const phaseModel& phase = this->fluid().mesh().template
    //     lookupObject<phaseModel>(this->thermo_->phasePropertyName("alpha"));

    //- Lable of water (H2O)
    const label specieI =
        this->thermo_->composition().species()[water];

            //- Get the specieStoichCoeff
    dimensionedScalar specieStoichCoeff
    (
        "stoichCoeff",
        dimless,
        rxnList_.found(water)
      ? rxnList_[water]/mag(rxnList_["e"])
      : 0.0
    );

    //- g/mol -> kg/mol
    const dimensionedScalar Wi
    (
        "W",
        dimMass/dimMoles,
        this->thermo_->composition().Wi(specieI)/1000
    );

    //- water production
    volScalarField wSpecie
    (
        i_*Wi*specieStoichCoeff/phaseModel::dimF
    );

    const word ionCatalystName = word(this->fluid().lookup("ionCatalystName"));
    label ionCatalystID = ionPhase.mesh().boundaryMesh().findPatchID(ionCatalystName);
    const polyPatch& ionCatalystPatch = ionPhase.mesh().boundaryMesh()[ionCatalystID];

    const word fluidCatalystName = word(this->fluid().lookup("fluidCatalystName"));
    label fluidCatalystID = fluidPhase.mesh().boundaryMesh().findPatchID(fluidCatalystName);
    const polyPatch& fluidCatalystPatch = fluidPhase.mesh().boundaryMesh()[fluidCatalystID];

    //- Dissolved water model
    surfaceDissolvedModel& dW = const_cast<surfaceDissolvedModel&>
    (
        ionPhase.mesh().template
        lookupObject<surfaceDissolvedModel>(surfaceDissolvedModel::modelName)
    );

    //- Water activity and water source/sink
    volScalarField& act = const_cast<volScalarField&>(dW.act());
    fvPatchField<scalar>& actPatch = act.boundaryField()[ionCatalystID];

    fvPatchField<scalar>& wSpeciePatch = wSpecie.boundaryField()[fluidCatalystID];

    //- Water mole fraction
    const volScalarField& XH2O = this->X(water);

    const volScalarField pSat = saturation_->pSat(this->thermo_->T());
    // Info << "Print pSat " << pSat << endl;

    const fvPatchField<scalar>& pSatPatch = pSat.boundaryField()[fluidCatalystID];


    //- Activity
    scalarField act0 =
    XH2O.boundaryField()[fluidCatalystID]
    * this->thermo_->p().boundaryField()[fluidCatalystID]
    / pSatPatch
    + 2.*(scalar(1) - this->boundaryField()[fluidCatalystID]);

    patchToPatchInterpolation ionToFluid
    (
        ionCatalystPatch,
        fluidCatalystPatch
    );

    patchToPatchInterpolation fluidToIon
    (
        fluidCatalystPatch,
        ionCatalystPatch
    );

    actPatch = fluidToIon.faceInterpolate(act0);

    //- Update the source/sink term dmdt
    dW.update(ionCatalystName);

    volScalarField& dmdt = const_cast<volScalarField&>(dW.dmdt());
    fvPatchField<scalar>& dmdtPatch = dmdt.boundaryField()[ionCatalystID];

    //- Relax
    // scalar iDmdtRelax(this->mesh().fieldRelaxationFactor("iDmdt")); 
    scalar iDmdtRelax = 1; 

    scalarField dmdtPatchOnFluid = ionToFluid.faceInterpolate(dmdtPatch);

    // TODO: Only consider single phase system
    //- Water is transferred between current phase and dissolved phase
    volScalarField& iDmdtWater = const_cast<volScalarField&>
        (this->iDmdt(water));

    fvPatchField<scalar>& iDmdtWaterPatch = iDmdtWater.boundaryField()[fluidCatalystID];

    // // Update Patch value for YiEqn
    // iDmdtWaterPatch = (1 - iDmdtRelax)*iDmdtWaterPatch
    // + iDmdtRelax*wSpeciePatch ;

    iDmdtWaterPatch = (1 - iDmdtRelax)*iDmdtWaterPatch
    + iDmdtRelax*(wSpeciePatch - dmdtPatchOnFluid*Wi.value());

    Info << "Print iDmdtWaterPatch for water " << average(iDmdtWaterPatch) << endl;
}

template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::SurfaceElectrochemicalReactingPhaseModel<BasePhaseModel>::R
(
    volScalarField& Y
) const
{
    FatalErrorInFunction
    << " R() function in surfaceElectroChemicalReaction is not used, rather than correctGasIdmdt "
    << exit(FatalError);
}

template<class BasePhaseModel>
void Foam::SurfaceElectrochemicalReactingPhaseModel<BasePhaseModel>::correctGasIdmdt
(
    const volScalarField& Y
) const
{

    word water("H2O");

    const label specieI =
        this->thermo_->composition().species()[Y.member()];

    //- kg/kmol -> kg/mol
    const dimensionedScalar Wi
    (
        "W",
        dimMass/dimMoles,
        this->thermo_->composition().Wi(specieI)/1000
    );

    dimensionedScalar specieStoichCoeff
    (
        "stoichCoeff",
        dimless,
        rxnList_.found(Y.member())
      ? rxnList_[Y.member()]/mag(rxnList_["e"])
      : 0.0
    );

    volScalarField wSpecie
    (
        i_*Wi*specieStoichCoeff/phaseModel::dimF
    );

    const word fluidCatalystName = word(this->fluid().lookup("fluidCatalystName"));
    label fluidCatalystID = this->mesh().boundaryMesh().findPatchID(fluidCatalystName);
    const polyPatch& fluidCatalystPatch = this->mesh().boundaryMesh()[fluidCatalystID];

    //- Relax
    // scalar iDmdtRelax(this->mesh().fieldRelaxationFactor("iDmdt")); 
    scalar iDmdtRelax(1); 

    volScalarField& iDmdtWater = const_cast<volScalarField&>(this->iDmdt(Y.member()));

    fvPatchField<scalar>& iDmdtWaterPatch = iDmdtWater.boundaryField()[fluidCatalystID];
    fvPatchField<scalar>& wSpeciePatch = wSpecie.boundaryField()[fluidCatalystID];

    iDmdtWaterPatch = (1 - iDmdtRelax)*iDmdtWaterPatch + iDmdtRelax*wSpeciePatch; 

    Info << "Print iDmdtWaterPatch for specie " << Y.member() << " " << average(iDmdtWaterPatch) << endl;
}

template<class BasePhaseModel>
void Foam::SurfaceElectrochemicalReactingPhaseModel<BasePhaseModel>::correctOverpotential()
{
    correctNernst();

    // ===== 1.) Calculate the overpotential from Butler Volmer equation on the electrode ====== //

    const regionType& fluidPhase = this->region
    (
        word(this->regions_.subDict("fluid").lookup("name"))
    );
    const regionType& electronPhase = this->region
    (
        word(this->regions_.subDict("electron").lookup("name"))
    );
    const regionType& ionPhase = this->region
    (
        word(this->regions_.subDict("ion").lookup("name"))
    );

    // Get the patch fields
    const word fluidCatalystName = word(this->fluid().lookup("fluidCatalystName"));
    label fluidCatalystID = fluidPhase.mesh().boundaryMesh().findPatchID(fluidCatalystName);
    const polyPatch& fluidCatalystPatch = fluidPhase.mesh().boundaryMesh()[fluidCatalystID];

    const word ionCatalystName = word(this->fluid().lookup("ionCatalystName"));
    label ionCatalystID = ionPhase.mesh().boundaryMesh().findPatchID(ionCatalystName);
    const polyPatch& ionCatalystPatch = ionPhase.mesh().boundaryMesh()[ionCatalystID];

    const word electronCatalystName = word(this->fluid().lookup("electronCatalystName"));
    label electronCatalystID = electronPhase.mesh().boundaryMesh().findPatchID(electronCatalystName);
    const polyPatch& electronCatalystPatch = electronPhase.mesh().boundaryMesh()[electronCatalystID];

    // patchToPatch interpolation
    patchToPatchInterpolation electronToFluid
    (
        electronCatalystPatch,
        fluidCatalystPatch
    );

    patchToPatchInterpolation ionToFluid
    (
        ionCatalystPatch,
        fluidCatalystPatch
    );
    patchToPatchInterpolation fluidToElectron
    (
        fluidCatalystPatch,
        electronCatalystPatch
    );

    patchToPatchInterpolation fluidToIon
    (
        fluidCatalystPatch,
        ionCatalystPatch
    );

    //- Access temperature and reactant
    const volScalarField& T = this->thermo_->T();
    const fvPatchField<scalar>& TPatch = T.boundaryField()[fluidCatalystID];

    // Get the alpha field
    const volScalarField& s = *this;
    const fvPatchField<scalar>& sPatch = s.boundaryField()[fluidCatalystID];

    //- Mixture mole fraction
    const volScalarField W(this->thermo_->W()/1000.0);
    const fvPatchField<scalar>& WPatch = W.boundaryField()[fluidCatalystID];
    
    //- Mixture density
    const volScalarField& rho= this->thermo_->rho();
    const fvPatchField<scalar>& rhoPatch = rho.boundaryField()[fluidCatalystID];

    //- Current exchange density
    fvPatchField<scalar>& i0Patch = this->i0_.boundaryField()[fluidCatalystID];

    //- Overpotential
    fvPatchField<scalar>& etaPatch = this->eta_.boundaryField()[fluidCatalystID];

    // Current density - from ionic region
    // Get the value of the current density, which was updated on the electronIon regions
    // and calculate the jump here on the fluid region 
    const volVectorField& iIonPhase = ionPhase.mesh().objectRegistry::lookupObject<volVectorField>
    (
        word(this->iNames_["ion"])
    );

    const fvPatchField<vector>& iIonPatch = iIonPhase.boundaryField()[ionCatalystID];
    scalarField iIonPatchSf = iIonPatch & (iIonPatch.patch().Sf() / iIonPatch.patch().magSf());
    scalarField iOnFluidPatch = mag(ionToFluid.faceInterpolate(iIonPatchSf));

    fvPatchField<scalar>& iFluidPatch = i_.boundaryField()[fluidCatalystID];

    Info << "Print average(iOnFluidPatch) " << average(iOnFluidPatch) << endl;
    iFluidPatch = iOnFluidPatch;

    Info << "Print average(iFluidPatch) " << average(iFluidPatch) << endl;

    //- Store values for all species
    scalarField coeff(TPatch.size(),1);

    //- Loop for all relative species
    forAllConstIter(dictionary, species_, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& dict0 = iter().dict();

            scalar ksi = readScalar(dict0.lookup("ksi"));
            scalar cRef = readScalar(dict0.lookup("cRef"));

            const volScalarField& X = this->X(name);
            const fvPatchField<scalar>& XPatch = X.boundaryField()[fluidCatalystID];

            coeff *= Foam::pow(XPatch*rhoPatch/WPatch/cRef, ksi);
        }
    }

    //- Loop for all relative liquid species
    forAllConstIter(dictionary, species_, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& dict0 = iter().dict();

            scalar ksi = readScalar(dict0.lookup("ksi"));
            scalar cRef = readScalar(dict0.lookup("cRef"));

            const volScalarField& X = this->X(name);
            const fvPatchField<scalar>& XPatch = X.boundaryField()[fluidCatalystID];

            coeff *= Foam::pow(XPatch*rhoPatch/WPatch/cRef, ksi);
        }
    }

    i0Patch = i0Value_.value()*coeff*Foam::pow(sPatch, this->gamma_);

    // Testing
    Info << "Print this->alpha_ " << this->alpha_ << endl;
    Info << "Print phaseModel::dimF.value() " << phaseModel::dimF.value() << endl;
    Info << "Print phaseModel::Rgas.value() " << phaseModel::Rgas.value() << endl;
    Info << "Print average(TPatch) " << average(TPatch) << endl;
    Info << "Print average(i0Patch) " << average(i0Patch) << endl;
    Info << "Print average(iOnFluidPatch) " << average(iOnFluidPatch) << endl;
    Info << "Print average(iFluidPatch) " << average(iFluidPatch) << endl;

    forAll(fluidCatalystPatch, faceI)
    {
            // newtonRaphson procedure

            //- electron transfer
            // //- cathode side, < 0
            // //- anode side, >0
            const scalar n = rxnList_["e"];
            const scalar sign = n/mag(n);

            A_ = n*this->alpha_*phaseModel::dimF.value()/phaseModel::Rgas.value()/TPatch[faceI];
            B_ = -n*(scalar(1) - this->alpha_)*phaseModel::dimF.value()/phaseModel::Rgas.value()/TPatch[faceI];

            //- Get eta with the root finder
            etaPatch[faceI] = this->eta(-1.2, 1.2, i0Patch[faceI], iOnFluidPatch[faceI]);
    }

    // Get the species
    const PtrList<volScalarField>& Y = this->Y();

    fvPatchField<scalar>& jumpPatch = this->jump_.boundaryField()[fluidCatalystID];

    //- Nernst field
    const fvPatchField<scalar>& nernstPatch = nernst_.boundaryField()[fluidCatalystID];

    jumpPatch = etaPatch + nernstPatch; // + etaContact; //+ etaConcentration ;

    Info << "Print etaPatch " << average(etaPatch) << endl;
    Info << "Print jumpPatch " << average(jumpPatch) << endl;

}


template<class BasePhaseModel>
void Foam::SurfaceElectrochemicalReactingPhaseModel<BasePhaseModel>::correctNernst()
{
    this->deltaH_ *= 0.0;
    this->deltaS_ *= 0.0;

    const word water = phaseModel::water;

    const volScalarField& p = this->thermo_->p();
    const volScalarField& T = this->thermo_->T();

    const word fluidCatalystName = word(this->fluid().lookup("fluidCatalystName"));
    label fluidCatalystID = this->mesh().boundaryMesh().findPatchID(fluidCatalystName);

    scalarField& nernstPatch(nernst_.boundaryField()[fluidCatalystID]);
    scalarField TCatalystPatch = T.boundaryField()[fluidCatalystID];
    scalarField pCatalystPatch = p.boundaryField()[fluidCatalystID];
    scalarField& deltaHPatch = this->deltaH_.boundaryField()[fluidCatalystID];
    scalarField& deltaSPatch = this->deltaS_.boundaryField()[fluidCatalystID];
    scalarField Qrxn(TCatalystPatch.size(), 1.0);
    scalarField deltaHSPatch(TCatalystPatch.size(), 0.0);
    scalarField pRef = pCatalystPatch/this->pRef_.value();

    forAllConstIter(HashTable<scalar>, rxnList_, iter)
    {
        const word& nameI = iter.key();

        if (nameI != "e") //&& nameI != water) TODO: For now only singlePhase electrochemistry
        {
            label speciesI = this->thermo_->composition().species()[nameI];
            const scalar Wi = this->thermo_->composition().Wi(speciesI)/1000.0;
            scalar stoiCoeffI = rxnList_[nameI];

            forAll(pCatalystPatch, faceI)
            {
                deltaHPatch[faceI] +=
                    stoiCoeffI*this->thermo_->composition().Ha(speciesI, pCatalystPatch[faceI], TCatalystPatch[faceI])*Wi;

                deltaSPatch[faceI] +=
                    stoiCoeffI*this->thermo_->composition().S(speciesI, pCatalystPatch[faceI], TCatalystPatch[faceI])*Wi;

                deltaHSPatch[faceI] +=
                    stoiCoeffI*this->thermo_->composition().S(speciesI, pCatalystPatch[faceI], TCatalystPatch[faceI])*Wi*TCatalystPatch[faceI];
            }

            const scalarField& X = this->X(nameI).boundaryField()[fluidCatalystID];

            Qrxn *= Foam::pow(Foam::max(X, 1.0e-6)*pRef, stoiCoeffI);
        }
    }

    nernstPatch = -(-(deltaHPatch - deltaHSPatch) - phaseModel::Rgas.value()*TCatalystPatch*Foam::log(Qrxn))/this->rxnList_["e"]/phaseModel::dimF.value();

    Info<< "Nernst " << fluidCatalystName
        << ": min = " << Foam::min(nernstPatch)
        << ", mean = " << Foam::average(nernstPatch)
        << ", max = " << Foam::max(nernstPatch)
        << endl;


}

template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::SurfaceElectrochemicalReactingPhaseModel<BasePhaseModel>::Qdot()
{
    tmp<volScalarField> q
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("q", this->name()),
                eta_.mesh().time().timeName(),
                eta_.mesh()
            ),
            eta_.mesh(),
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0)
        )
    );

    // if (nernst_ == nullptr)
    // {
    //     return q;
    // }

        //- get the sub models
    const regionType& fluidPhase = region
    (
        word(regions_.subDict("fluid").lookup("name"))
    );

    const word fluidCatalystName = word(this->fluid().lookup("fluidCatalystName"));
    label fluidCatalystID = fluidPhase.mesh().boundaryMesh().findPatchID(fluidCatalystName);
    const polyPatch& fluidCatalystPatch = fluidPhase.mesh().boundaryMesh()[fluidCatalystID];

    const volScalarField& i = i_;
    const fvPatchField<scalar>& iPatch = i.boundaryField()[fluidCatalystID];

    const volScalarField& T = this->thermo().T();
    const fvPatchField<scalar>& TPatch = T.boundaryField()[fluidCatalystID];

    const volScalarField& eta = eta_;
    const fvPatchField<scalar>& etaPatch = eta.boundaryField()[fluidCatalystID];

    const volScalarField& deltaS = this->deltaS();
    const fvPatchField<scalar>& deltaSPatch = deltaS.boundaryField()[fluidCatalystID];

    scalar n = rxnList_["e"];
    scalar sign = n/mag(n);

    tmp<scalarField> qPatch = iPatch * (sign*etaPatch - TPatch*deltaSPatch
                              /mag(rxnList_["e"])/phaseModel::dimF.value());

    // const labelList& fluidCatalystFaceCells = fluidCatalystPatch.faceCells(); 

    forAll(fluidCatalystPatch, faceI)
    {
        label faceCelli = fluidCatalystPatch.faceCells()[faceI];
        q()[faceCelli] = qPatch()[faceI]* (mag(fluidCatalystPatch.faceAreas()[faceI]) / this->mesh().V()[faceCelli]); ;
    }

    Qdot_ = q();

    // Info << "Heat flux due to the electrochemical reactions " << average(q()) << endl;                     

    return q;
}


template<class BasePhaseModel>
Foam::scalar
Foam::SurfaceElectrochemicalReactingPhaseModel<BasePhaseModel>::eta
(
    scalar etaMin,
    scalar etaMax,
    scalar i0,
    scalar i
) const
{
    scalar eps(GREAT);
    scalar eta0(0);

    //- Newton-Raphson method
    do
    {
        scalar etaOld = eta0;

        scalar g =
            i
          - i0
            *  (
                   Foam::exp(A_*eta0)
                 - Foam::exp(B_*eta0)
               );

        scalar f =
          - i0
            *  (
                   A_*Foam::exp(A_*eta0)
                 - B_*Foam::exp(B_*eta0)
               );

        eta0 -= g/f;

        eta0 = Foam::max(eta0, etaMin);
        eta0 = Foam::min(eta0, etaMax);

        eps = Foam::min(eps, Foam::mag(etaOld - eta0));

    } while (eps > 1.0e-10);

    return eta0;
}
// ************************************************************************* //
