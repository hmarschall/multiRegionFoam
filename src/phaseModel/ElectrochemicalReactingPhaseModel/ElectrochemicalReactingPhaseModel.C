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

#include "ElectrochemicalReactingPhaseModel.H"
// #include "phaseSystem.H"
#include "dissolvedModel.H"


const Foam::dimensionedScalar Rgas("Rgas", Foam::dimensionSet(1,2,-2,-1,-1,0,0), 8.314462); // [J/mol/K]
const Foam::dimensionedScalar dimF("dimF", Foam::dimensionSet(0,0,1,0,-1,1,0), 96485.33212); // [C/mol]

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::ElectrochemicalReactingPhaseModel<BasePhaseModel>::ElectrochemicalReactingPhaseModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    BasePhaseModel(dict, mesh),
    eta_
    (
        IOobject
        (
            "eta",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar
        (
            "eta",
            dimensionSet(1, 2, -3, 0, 0, -1, 0),
            0.0
        )
    ),
    zoneName_(dict.lookup("zoneName")),
    regions_(dict.subDict("regions")),
    saturation_(saturationModel::New(dict.subDict("saturation"), this->mesh())),
    phiNames_(dict.lookup("phiNames")),
    relax_(dict.lookupOrDefault<scalar>("relax", 1.0)),
    alpha_(dict.lookupOrDefault<scalar>("alpha", 0.5)),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 1.0)),
    j0Value_("j0", dimCurrent/dimVolume, readScalar(dict.lookup("j0"))),
    j_
    (
        IOobject
        (
            "j",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar
        (
            "j",
            dimCurrent/dimVolume,
            0.0
        ),
        zeroGradientFvPatchScalarField::typeName
    ),
    j0_
    (
        IOobject
        (
            "j0",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar
        (
            "j0",
            dimCurrent/dimVolume,
            0.0
        ),
        zeroGradientFvPatchScalarField::typeName
    ),
    EA_("EA", dimEnergy/dimMoles, readScalar(dict.lookup("EA"))),
    Tref_("Tref", dimTemperature, readScalar(dict.lookup("Tref"))),
    rxnList_(dict.lookup("RxnList")),
    // Nernst variables
    nernst_
    (
        IOobject
        (
            "nernst",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
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
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
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
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar
        (
            "deltaS",
            dimEnergy/dimMoles/dimTemperature,
            0.0
        )
    ),
    residualY_(dict.lookupOrDefault<scalar>("residualY", 1.0e-6)),
    pRef_("pRef", dimPressure, readScalar(dict.lookup("pRef"))),
    species_(dict.subDict("species"))
{
    // Info << "Breakpoint ElectrochemicalReactingPhaseModel.C " << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::ElectrochemicalReactingPhaseModel<BasePhaseModel>::~ElectrochemicalReactingPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::ElectrochemicalReactingPhaseModel<BasePhaseModel>::correct()
{
    //- Activation overpotential model update
    correctOverpotential();

    //- Get sub regions
    //- Refer to regionType
    //- Including: fluid, electron (BPP + GDL + CL), ion (CLs + membrane)
    const regionType& fluidPhase = this->mesh().time().template
                 lookupObject<regionType>(word(this->dict().subDict("fluid").lookup("name")));

    const regionType& ionPhase = this->mesh().time().template 
                lookupObject<regionType>(word(this->dict().subDict("ion").lookup("name")));

    // // Get the phase System
    //  const phaseSystem& phaseSys = this->mesh().template
    //  lookupObject<phaseSystem>(phaseSystem::propertiesName);

    // //- Get the present phase Model from fluid phase
    // const phaseModel& phase = this->mesh().template
    //     lookupObject<phaseModel>(this->thermo_->phasePropertyName("alpha"));

    word water("H2O");

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
        this->thermo_->composition().W(specieI)/1000
    );

    //- water production
    volScalarField wSpecie
    (
        j_*Wi*specieStoichCoeff/dimF
    );

    //- Dissolved water model
    dissolvedModel& dW = const_cast<dissolvedModel&>
    (
        ionPhase.template
        lookupObject<dissolvedModel>(dissolvedModel::modelName)
    );

    //- Water activity and water source/sink
    scalarField& act = const_cast<volScalarField&>(dW.act());
    scalarField& dmdt = const_cast<volScalarField&>(dW.dmdt());

    //- Water mole fraction
    const scalarField& xH2O = this->X(water);

    //- Activity
    scalarField act0 =
        xH2O
      * this->thermo_->p()
      / saturation_->pSat(this->thermo_->T())()
      + 2*(scalar(1) - this->internalField());

    //- Only consider catalyst zone
    label znId = fluidPhase.cellZones().findZoneID(zoneName_);
    const labelList& cells = fluidPhase.cellZones()[znId];

    //- Update water activity
    forAll(cells, cellI)
    {
        //- get cell IDs
        label fluidId = cells[cellI];
        label ionId = ionPhase.cellMap()[fluidPhase.cellMapIO()[fluidId]];

        //- Water activity: act = p(H2O)/pSat + 2*sat
        act[ionId] = act0[fluidId];
    }

    //- Update the source/sink term dmdt
    dW.update(zoneName_);

    //- Relax
    // scalar iDmdtRelax(this->mesh().fieldRelaxationFactor("iDmdt")); 
    scalar iDmdtRelax = 1; 

    //- Update the water production in phase model
    forAll(cells, cellI)
    {
        //- get cell IDs
        label fluidId = cells[cellI];
        label ionId = ionPhase.cellMap()[fluidPhase.cellMapIO()[fluidId]];

        if (dissolved_)
        {
            dmdt[ionId] += wSpecie[fluidId]/Wi.value();
        }

        // if (phaseSys.isSinglePhase() || !eta_->phaseChange())
        // {
            //- Water is transferred between current phase and dissolved phase
            volScalarField& iDmdtWater = const_cast<volScalarField&>
                (this->iDmdt(water));

            iDmdtWater[fluidId] = (1 - iDmdtRelax)*iDmdtWater[fluidId]
                + iDmdtRelax*(wSpecie[fluidId] - dmdt[ionId]*Wi.value());
        // }
        // else
        // {
        //     //- Get the name of the other phase
        //     const word name1 = Pair<word>
        //     (
        //         phaseSys.phases()[0].name(),
        //         phaseSys.phases()[1].name()
        //     ).other(this->name());

        //     //- Water is transferred between the other phase and dissolved phase
        //     scalarField& iDmdtWater = const_cast<volScalarField&>
        //         (phaseSys.phases()[name1].iDmdt(water));

        //     iDmdtWater[fluidId] = (1 - iDmdtRelax)*iDmdtWater[fluidId]
        //         + iDmdtRelax*(wSpecie[fluidId] - dmdt[ionId]*Wi.value());
        // }
    }
}

template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::ElectrochemicalReactingPhaseModel<BasePhaseModel>::R
(
    volScalarField& Y
) const
{

    word water("H2O");

    if (Y.member() == water)
    {
        // Jap, ist auch immer 0
        // Info << "Print gaseous water idmdt (müsste immer 0 sein) " << this->iDmdt(water) << endl;
        return this->iDmdt(water) + fvm::Sp(this->iDmdt(water)*0.0, Y);
        // return -fvm::SuSp(-this->iDmdt(water),Y);
    }

    const label specieI =
        this->thermo_->composition().species()[Y.member()];

    //- kg/kmol -> kg/mol
    const dimensionedScalar Wi
    (
        "W",
        dimMass/dimMoles,
        this->thermo_->composition().W(specieI)/1000
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
        j_*Wi*specieStoichCoeff/dimF
    );


    //- Relax
    // scalar iDmdtRelax(this->mesh().fieldRelaxationFactor("iDmdt")); 
    scalar iDmdtRelax = 1;
// vsf iDmdt in multiComponentFickPhaseModel ist ja nur, um dieses Feld als VTK dann ausschreiben und anshenen zu können
    volScalarField& iDmdtWater = const_cast<volScalarField&>(this->iDmdt(Y.member()));      

    iDmdtWater = (1 - iDmdtRelax)*iDmdtWater + iDmdtRelax*wSpecie; 

    return iDmdtWater + fvm::Sp(wSpecie*0.0, Y);
    // return -fvm::SuSp(-this->iDmdt(water),Y);

}

template<class BasePhaseModel>
void Foam::ElectrochemicalReactingPhaseModel<BasePhaseModel>::correctOverpotential()
{
    correctNernst();

    //- Get sub regions
    //- Refer to regionType
    //- Including: fluid, electron (BPP + GDL + CL), ion (CLs + membrane)
    // const regionType& fluidPhase = this->region
    // (
    //     word(this->regions_.subDict("fluid").lookup("name"))
    // );
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

    //- Source/sink terms for electron and ion fields
    //- See names in regions/electronIon/electronIon.C
    scalarField& SE = const_cast<volScalarField&>
    (
        electronPhase.template lookupObject<volScalarField>("J")
    );
    scalarField& SI = const_cast<volScalarField&>
    (
        ionPhase.template lookupObject<volScalarField>("J")
    );

    //- Potential fields
    const scalarField& phiE = electronPhase.template
        lookupObject<volScalarField>
        (
            word(this->phiNames_["electron"])
        );
    const scalarField& phiI = ionPhase.template
        lookupObject<volScalarField>
        (
            word(this->phiNames_["ion"])
        );

    //- Reference
    scalarField& eta = eta_;
    //- Nernst field
    scalarField& nernst = nernst_;
    //- Current density
    scalarField& j = j_;

    //- Access temperature and reactant
    const scalarField& T = this->thermo_->T();
    const scalarField& s = *this;

    //- Mixture mole fraction
    const scalarField W(this->thermo_->W()/1000.0);
    //- Mixture density
    const scalarField& rho= this->thermo_->rho();
    //- Store values for all species

    scalarField coeff(this->mesh().nCells(), 1.0);
    //- Loop for all relative species
    forAllConstIter(dictionary, this->species_, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& dict0 = iter().dict();

            scalar ksi = readScalar(dict0.lookup("ksi"));
            scalar cRef = readScalar(dict0.lookup("cRef"));

            const scalarField& X = this->X(name);

            coeff *= Foam::pow(X*rho/W/cRef, ksi);
        }
    }

    //- Only consider catalyst zone
    label znId = fluidPhase.cellZones().findZoneID(this->zoneName_);
    const labelList& cells = fluidPhase.cellZones()[znId];

    //- electron transfer
    //- cathode side, < 0
    //- anode side, >0
    scalar n = rxnList_["e"];
    scalar sign = n/mag(n);

    //- The total current: volume averaged
    scalar Rj(0.0);

    volScalarField& j0 = j0_;

    forAll(cells, cellI)
    {
        //- get cell IDs
        label fluidId = cells[cellI];
        label electronId = electronPhase.cellMap()[fluidPhase.cellMapIO()[fluidId]];
        label ionId = ionPhase.cellMap()[fluidPhase.cellMapIO()[fluidId]];

        //- activation overpotential
        eta[fluidId] = 
        (
            eta[fluidId]*(scalar(1) - relax_)
          + (phiE[electronId] - phiI[ionId] - nernst[fluidId])
          * relax_
        );


        j0[fluidId] = j0Value_.value()
        // *Foam::pow(s[fluidId], this->gamma_)
        // *Foam::pow(1-s[fluidId], this->gamma_)  // For the case of an electrolyzer, i0 decreases if more gas is around the electrode surface
        *Foam::exp(EA_.value()/Rgas.value()*((1/Tref_.value()) - (1/T[fluidId])))
        *coeff[fluidId]; //*Foam::pow(s[fluidId], this->gamma_)* 

        // Info << "Print eta[fluidId] " << eta[fluidId] << endl;

        //- Buttler-volmer relation
        j[fluidId] = Foam::max
        (
            j0[fluidId]*
            (
                Foam::exp(n*alpha_*dimF.value()*eta[fluidId]/Rgas.value()/T[fluidId])
              - Foam::exp(-n*(scalar(1) - alpha_)*dimF.value()*eta[fluidId]/Rgas.value()/T[fluidId])
            )
            ,
            scalar(0)
        );

        //- anode side: SE = Rj, SI = -Rj
        //- cathode side: SE = -Rj, SI = Rj
        SE[electronId] = -sign*j[fluidId];
        SI[ionId] = -SE[electronId];

        Rj += fluidPhase.V()[fluidId] * j[fluidId];
    }

    reduce(Rj, sumOp<scalar>());

    Info << "Total current (A) at " << this->zoneName_
         << ": " << Rj << endl;
}


template<class BasePhaseModel>
void Foam::ElectrochemicalReactingPhaseModel<BasePhaseModel>::correctNernst()
{

    this->deltaH_ *= 0.0;
    this->deltaS_ *= 0.0;

    const word water = "H2O";

    scalarField Qrxn(this->size(), 1.0);
    scalarField deltaHS(this->size(), 0.0);

    const scalarField& p = this->thermo_->p();
    const scalarField& T = this->thermo_->T();

    scalarField pRef = p/this->pRef_.value();

    forAllConstIter(HashTable<scalar>, this->rxnList_, iter)
    {
        const word& nameI = iter.key();

        if (nameI != "e") //&& nameI != water) TODO: For now only singlePhase electrochemistry
        {
            label speciesI = this->thermo_->composition().species()[nameI];
            const scalar Wi = this->thermo_->composition().W(speciesI)/1000.0;
            scalar stoiCoeffI = rxnList_[nameI];

            forAll(this->deltaH_, cellI)
            {
                this->deltaH_[cellI] +=
                    stoiCoeffI*this->thermo_->composition().Ha(speciesI, p[cellI], T[cellI])*Wi;

                this->deltaS_[cellI] +=
                    stoiCoeffI*this->thermo_->composition().S(speciesI, p[cellI], T[cellI])*Wi;

                deltaHS[cellI] +=
                    stoiCoeffI*this->thermo_->composition().S(speciesI, p[cellI], T[cellI])*Wi*T[cellI];
            }

            const scalarField& X = this->X(nameI);

            Qrxn *= Foam::pow(Foam::max(X, 1.0e-6)*pRef, stoiCoeffI);
        }
    }

    // nernst_ = -(-(this->deltaH_ - deltaHS) - Rgas*T*Foam::log(Qrxn))/this->rxnList_["e"]/F;

    Info<< "Nernst " << this->mesh().name()
        << ": min = " << Foam::min(this->internalField())
        << ", mean = " << Foam::average(this->internalField())
        << ", max = " << Foam::max(this->internalField())
        << endl;


}
// ************************************************************************* //
