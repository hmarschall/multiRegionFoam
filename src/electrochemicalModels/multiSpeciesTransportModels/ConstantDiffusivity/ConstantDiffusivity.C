/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "ConstantDiffusivity.H"
#include "fvmLaplacian.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

// template<class Thermo>
// void Foam::multiSpeciesTransportModels::ConstantDiffusivity<Thermo>::updateCoefficients()
// {

//     // const volScalarField& rho = this->phase_.thermo().rho();
//     const volScalarField& T = this->thermo_.T();
//     // const dimensionedScalar Tref("Tref", dimTemperature, this->getScalar("Tref"));

//     this->DijModel_().update();

//     // forAll(this->species_, i)
//     // {
//     //     forAll(this->species_, j)
//     //     {
//     //         Info << "Print this->species_[i] and this->species_[j] " << this->species_[i] << "and " << this->species_[j] << endl;
//     //         D_[i] = rho*this->Dij(i,j);
//     //         Info << "Print gAverage(D_[i]) " << gAverage(D_[i]) << endl;
//     //     }
//     // }
//     forAll(this->species_, i)
//     {
//             // Info << "Print this->species_[i] and this->species_[j] " << this->species_[i] << endl;
//             const dimensionedScalar EA("EA", dimEnergy/dimMoles, this->activationEnergyList_[this->species_[i]]);
//             const dimensionedScalar Tref("Tref", dimTemperature, this->TrefList_[this->species_[i]]);

//             D_[i] = this->Dij(i,i) *Foam::exp(EA/Rgas_*((1/Tref) - (1/T)));

//             // const volScalarField arh = this->Dij(i,i) *Foam::exp(EA/Rgas_*((1/Tref) - (1/T))) ;
//             // Info << "Print gAverage(this->Dij(i,i)) " << gAverage(this->Dij(i,i)) << endl;
//             // Info << "Print gAverage(this->Dij(i,i) *Foam::exp(EA/Rgas_*((1/Tref) - (1/T)))) " << gAverage(arh) << endl;
//             // Info << "Print rho " << rho << endl;
//             // Info << "Print Tref " << Tref << endl;
//             // Info << "Print gAverage(rho) " << gAverage(rho) << endl;
//             // Info << "Print i " << i << endl;
//             // Info << "Print gAverage(D_[i]) " << gAverage(D_[i]) << endl;
//     }

//     label speciesIndexOH = this->thermo_.composition().species()["OH-"];
//     label speciesIndexK = this->thermo_.composition().species()["K+"];
//     volScalarField& tOH = t_[speciesIndexOH];
//     volScalarField& tK = t_[speciesIndexK];
//     volScalarField& DOH = D_[speciesIndexOH]; 
//     volScalarField& DK = D_[speciesIndexK];
//     scalar zOH = this->valenceList_["OH-"];
//     scalar zK = this->valenceList_["K+"];

//     tOH = zOH*DOH /(zOH*DOH - zK*DK); 
//     tK = zK*DK /(zK*DK - zOH*DOH); 

// } 


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::multiSpeciesTransportModels::ConstantDiffusivity<Thermo>::ConstantDiffusivity
(
    const phaseModel& phase
)
:
    MultiSpeciesTransportModel<Thermo>(phase),
    X_(phase.X()),
    // constantDiffusivityDict_(this->subDict(typeName + "Coeff"))
    constantDiffusivityDict_(this->subDict("constantDiffusionCoeff")),
    activationEnergyList_
    (
        constantDiffusivityDict_.lookup("activationEnergyList")
    ),
    TrefList_
    (
        constantDiffusivityDict_.lookup("TrefList")
    ),
    Dm_(this->thermo_.composition().species().size()),
    DRefm_(this->thermo_.composition().species().size()),
    eps_
    (
        IOobject
        (
            "porosityFieldConstDiff",
            phase.mesh().time().timeName(),
            phase.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase.mesh(),
        dimensionedScalar("epsDiff", dimless, 1)
    ),
    tau_
    (
        IOobject
        (
            "tortuosityFieldConstDiff",
            phase.mesh().time().timeName(),
            phase.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase.mesh(),
        dimensionedScalar("tauDiff", dimless, 1)
    ),
    useArrh_(constantDiffusivityDict_.lookupOrDefault<Switch>("useArhenius", false)),
    Rgas_(dimensionedScalar("Rgas", dimensionSet(1,2,-2,-1,-1,0,0), 8.314462))
{
    const PtrList<volScalarField>& Y = phase.Y();

    Dm_.setSize(Y.size());
    // t_.setSize(Y.size());
    
    forAll(Dm_, i)
    {
        DRefm_.set
        (
            i, 
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("DRefm", Y[i].member()),
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh_,
                dimensionedScalar("DRefm", dimensionSet(0,2,-1,0,0,0,0), 
                                    readScalar(constantDiffusivityDict_.lookup(Y[i].member())))
            )
        );

        Dm_.set
        (
            i, 
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("Dm", Y[i].member() + "." + phase.name()),
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh_,
                dimensionedScalar("Dm", dimensionSet(0,2,-1,0,0,0,0), SMALL)
            )
        );
    }

    // Get porosity, tortuosity and pore diameter from porousZones
    forAll(this->pZones_, zoneI)
    {
        const word zoneName = this->pZones_[zoneI].name();
        label znId = this->mesh_.cellZones().findZoneID(zoneName);
        const labelList& cells = this->mesh_.cellZones()[znId];

        forAll(cells, cellI)
        {
            label porousID = cells[cellI];

            eps_[porousID] = this->pZones_.operator[](zoneI).porosity();
            tau_[porousID] = this->pZones_.operator[](zoneI).tau();
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::multiSpeciesTransportModels::ConstantDiffusivity<Thermo>::~ConstantDiffusivity()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


template<class Thermo>
Foam::tmp<Foam::volScalarField> Foam::multiSpeciesTransportModels::ConstantDiffusivity<Thermo>::DEff
(
    const volScalarField& Yi
) const
{
    const basicSpecieMixture& composition = this->thermo_.composition();

    Info << "Print Yi.member() " << Yi.member() << endl;
    Info << "Print Dm_[composition.index(Yi)] " << gAverage(Dm_[composition.index(Yi)]) << endl;

    return this->thermo_.rho()*Dm_[composition.index(Yi)];
}


template<class Thermo>
Foam::tmp<Foam::scalarField> Foam::multiSpeciesTransportModels::ConstantDiffusivity<Thermo>::DEff
(
    const volScalarField& Yi,
    const label patchi
) const
{
    const basicSpecieMixture& composition = this->thermo_.composition();

    return
        this->thermo_.rho()().boundaryField()[patchi]
       *Dm_[composition.index(Yi)].boundaryField()[patchi];
}


template<class Thermo>
Foam::tmp<Foam::surfaceScalarField> Foam::multiSpeciesTransportModels::ConstantDiffusivity<Thermo>::q() const
{

    // tmp<surfaceScalarField> tmpq
    // (
    //     surfaceScalarField::New
    //     (
    //         "q",
    //        -fvc::interpolate(this->phase_*this->kappaEff())
    //        *fvc::snGrad(this->thermo_.T())
    //     )
    // );

    tmp<surfaceScalarField> tmpq
    (
        new surfaceScalarField
        (
            IOobject
            (
                "q",
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            -fvc::interpolate(this->phase_*this->kappaEff())
            *fvc::snGrad(this->thermo_.T())
        )
    );

    const basicSpecieMixture& composition = this->thermo_.composition();
    const PtrList<volScalarField>& Y = composition.Y();

    if (Y.size())
    {
        tmp<surfaceScalarField> tsumJ
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "sumJ",
                    this->mesh_.time().timeName(),
                    this->mesh_
                ),
                this->mesh_,
                dimensionedScalar("sumJ", dimMass/dimArea/dimTime, 0)
            )
        );

        surfaceScalarField& sumJ = tsumJ();

        tmp<surfaceScalarField> tsumJh
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "sumJh",
                    this->mesh_.time().timeName(),
                    this->mesh_
                ),
                this->mesh_,
                dimensionedScalar("sumJh", dimMass/dimArea/dimTime*dimForce*dimLength/dimMass, 0)
            )
        );

        surfaceScalarField& sumJh = tsumJh();

        forAll(Y, i)
        {
            if (i != composition.defaultSpecie())
            {
                const volScalarField hi
                (
                    composition.Hs(i, this->thermo_.p(), this->thermo_.T())
                );

                const surfaceScalarField ji(this->j(Y[i]));
                sumJ += ji;

                sumJh += ji*fvc::interpolate(hi);
            }
        }

        {
            const label i = composition.defaultSpecie();

            const volScalarField hi
            (
                composition.Hs(i, this->thermo_.p(), this->thermo_.T())
            );

            sumJh -= sumJ*fvc::interpolate(hi);
        }

        tmpq() += sumJh;
    }

    return tmpq;
}


template<class Thermo>
Foam::tmp<Foam::fvScalarMatrix> Foam::multiSpeciesTransportModels::ConstantDiffusivity<Thermo>::divq
(
    volScalarField& he
) const
{
    tmp<fvScalarMatrix> tmpDivq
    (
        fvm::Su
        (
            -fvc::laplacian(this->phase_*this->kappaEff(), this->thermo_.T()),
            he
        )
    );

    const basicSpecieMixture& composition = this->thermo_.composition();
    const PtrList<volScalarField>& Y = composition.Y();

    tmpDivq() -=
        correction(fvm::laplacian(this->phase_*this->alphaEff(), he));

        tmp<surfaceScalarField> tsumJ
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "sumJ",
                    this->mesh_.time().timeName(),
                    this->mesh_
                ),
                this->mesh_,
                dimensionedScalar("sumJ", dimMass/dimArea/dimTime, 0)
            )
        );

        surfaceScalarField& sumJ = tsumJ();

        tmp<surfaceScalarField> tsumJh
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "sumJh",
                    this->mesh_.time().timeName(),
                    this->mesh_
                ),
                this->mesh_,
                dimensionedScalar("sumJh", dimMass/dimArea/dimTime*dimForce*dimLength/dimMass, 0)
            )
        );

        surfaceScalarField& sumJh = tsumJh();

    forAll(Y, i)
    {
        if (i != composition.defaultSpecie())
        {
            const volScalarField hi
            (
                composition.Hs(i, this->thermo_.p(), this->thermo_.T())
            );

            const surfaceScalarField ji(this->j(Y[i]));
            sumJ += ji;

            sumJh += ji*fvc::interpolate(hi);
        }
    }

    {
        const label i = composition.defaultSpecie();

        const volScalarField hi
        (
            composition.Hs(i, this->thermo_.p(), this->thermo_.T())
        );

        sumJh -= sumJ*fvc::interpolate(hi);
    }

    tmpDivq() += fvc::div(sumJh*he.mesh().magSf());

    return tmpDivq;
}


template<class Thermo>
Foam::tmp<Foam::surfaceScalarField> Foam::multiSpeciesTransportModels::ConstantDiffusivity<Thermo>::j
(
    const volScalarField& Yi
) const
{
    tmp<surfaceScalarField> tj
    (
        new surfaceScalarField
        (
            IOobject
            (
                 "j" + Yi.name(),
                 this->mesh_.time().timeName(),
                 this->mesh_
            ),
            this->mesh_,
            dimensionedScalar("j(" + Yi.name() + ')', dimensionSet(1,-2,-1,0,0,0,0), 0)
        )
    );

    tj = -fvc::interpolate(this->phase_*DEff(Yi))*fvc::snGrad(Yi); 

    return tj;
}


template<class Thermo>
Foam::tmp<Foam::fvScalarMatrix> Foam::multiSpeciesTransportModels::ConstantDiffusivity<Thermo>::divj
(
    volScalarField& Yi
) const
{
    return -fvm::laplacian(this->phase_*DEff(Yi), Yi);
}


template<class Thermo>
void Foam::multiSpeciesTransportModels::ConstantDiffusivity<Thermo>::correct()
{
    //- Consider temperature dependency for the diffusion coefficient

    const basicSpecieMixture& composition = this->thermo_.composition();
    const PtrList<volScalarField>& Y = composition.Y();
    const volScalarField& T = this->thermo_.T();

    forAll(Y, i)
    {
        const dimensionedScalar EA("EA", dimEnergy/dimMoles, this->activationEnergyList_[Y[i].member()]);
        const dimensionedScalar Tref("Tref", dimTemperature, this->TrefList_[Y[i].member()]);

        // Including:
        // - Impact of alpha != 1 through phase
        // - Impact of porous media
        // - Temperature dependency
        if(useArrh_)
        {
            // Dm_[i] =  this->phase_*(eps_/tau_)*DRefm_[i]*Foam::exp(EA/Rgas_*((1/Tref) - (1/T)));
            Dm_[i] =  (eps_/tau_)*DRefm_[i]*Foam::exp(EA/Rgas_*((1/Tref) - (1/T)));
            // Info << "Enter Arrhenius constDiffusivity " << endl;
        }
        else
        {
            // Dm_[i] =  this->phase_*(eps_/tau_)*DRefm_[i];
            Dm_[i] =  (eps_/tau_)*DRefm_[i];
            // Info << "Enter Normal constDiffusivity " << endl;
        }
    }
}

template<class Thermo>
void Foam::multiSpeciesTransportModels::ConstantDiffusivity<Thermo>::updateBinaryDiffusionCoeff()
{
    //- Dummy
}

template<class Thermo>
const Foam::volScalarField& Foam::multiSpeciesTransportModels::ConstantDiffusivity<Thermo>::Dij(label i, label j) const
{
    //- Dummy

    tmp<volScalarField> Dij
    (
        new volScalarField
        (
            IOobject
            (
                "Dij", this->phase_.name(),
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar("Dij", dimensionSet(1,-2,-1,0,0,0,0), 0)
        )
    );

    return Dij;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
