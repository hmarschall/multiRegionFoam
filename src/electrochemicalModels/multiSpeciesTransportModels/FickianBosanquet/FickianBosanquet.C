/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    along with OpenFOAM.  If not, see <http://www.gnthis->org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "FickianBosanquet.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvmSup.H"
#include "surfaceInterpolate.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::multiSpeciesTransportModels::FickianBosanquet<Thermo>::FickianBosanquet
(
    const phaseModel& phase
)
:
    MultiSpeciesTransportModel<Thermo>(phase),
    dp_(this->pZones_.size()),
    eps_(this->pZones_.size()),
    tau_(this->pZones_.size()),
    Dm_(this->thermo_.composition().species().size()),
    DmK_(this->thermo_.composition().species().size()),
    DK_(this->thermo_.composition().species().size()),
    X_(phase.X()),
    DijModels_(0.5*this->thermo_.composition().species().size()*(this->thermo_.composition().species().size()+1)),
    Dij_(DijModels_.size())
{
    const basicSpecieMixture& composition = this->thermo_.composition();
    const PtrList<volScalarField>& Y = composition.Y();

    forAll(Y, i)
    {
        DK_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "DK_" + composition.Y()[i].name(),
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedScalar("DK", dimensionSet(0, 2, -1, 0, 0), 0)
            )
        );

        // Write DMixture as a field
        Dm_.set
        (
            i, 
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("Dm", Y[i].member()),
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh_,
                dimensionedScalar("Dm", dimensionSet(0,2,-1,0,0,0,0), Foam::SMALL)
            )
        );

        // Write DMixture as a field
        DmK_.set
        (
            i, 
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("DmK", Y[i].member()),
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh_,
                dimensionedScalar("DmK", dimensionSet(0,2,-1,0,0,0,0), Foam::SMALL)
            )
        );


        for(label j=i; j < Y.size(); j++)
        {
            label k = Y.size()*i+j-0.5*i*(i+1);
            
            DijModels_.set
            (
                k,
                binaryDiffusivityModel<Thermo>::New
                (
                    this->thermo_,
                    Y[i].member(),
                    Y[j].member(),
                    *this,
                    this->phase_
                )
            );

            Dij_.set
            (
                k,
                new volScalarField
                (
                   DijModels_[k].D() 
                )
            );  
        }
    }

    // Get porosity, tortuosity and pore diameter from porousZones
    forAll(this->pZones_, i)
    {
        Info << "Print i " << i << endl; 
        eps_[i] = this->pZones_.operator[](i).porosity();
        tau_[i] = this->pZones_.operator[](i).tau();
        dp_[i] = this->pZones_.operator[](i).poreDiametre();

        Info << "Print eps_[i] " << eps_[i] << endl;
        Info << "Print tau_[i] " << tau_[i] << endl;
        Info << "Print dp_[i] " << dp_[i] << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::multiSpeciesTransportModels::FickianBosanquet<Thermo>::~FickianBosanquet()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Thermo>
Foam::tmp<Foam::volScalarField> Foam::multiSpeciesTransportModels::FickianBosanquet<Thermo>::DEff
(
    const volScalarField& Yi
) const
{
    const basicSpecieMixture& composition = this->thermo_.composition();

    return this->thermo_.rho()*DmK_[composition.index(Yi)];
}


template<class Thermo>
Foam::tmp<Foam::scalarField> Foam::multiSpeciesTransportModels::FickianBosanquet<Thermo>::DEff
(
    const volScalarField& Yi,
    const label patchi
) const
{
    const basicSpecieMixture& composition = this->thermo_.composition();

    return 
        this->thermo_.rho()().boundaryField()[patchi]
        *DmK_[composition.index(Yi)].boundaryField()[patchi];
}


template<class Thermo>
Foam::tmp<Foam::surfaceScalarField> Foam::multiSpeciesTransportModels::FickianBosanquet<Thermo>::q() const
{


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
Foam::tmp<Foam::fvScalarMatrix> Foam::multiSpeciesTransportModels::FickianBosanquet<Thermo>::divq
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
Foam::tmp<Foam::surfaceScalarField> Foam::multiSpeciesTransportModels::FickianBosanquet<Thermo>::j
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

    tj = -fvc::interpolate(this->phase_*this->thermo_.rho()*DEff(Yi))*fvc::snGrad(Yi);      // This would call Fickian:.DEff (if part of template) 

    return tj;
}


template<class Thermo>
Foam::tmp<Foam::fvScalarMatrix> Foam::multiSpeciesTransportModels::FickianBosanquet<Thermo>::divj
(
    volScalarField& Yi
) const
{
    return -fvm::laplacian(this->phase_*DEff(Yi), Yi);
}


template<class Thermo>
void Foam::multiSpeciesTransportModels::FickianBosanquet<Thermo>::correct()
{

    updateBinaryDiffusionCoeff();

    const basicSpecieMixture& composition = this->thermo_.composition();
    const PtrList<volScalarField>& Y = composition.Y();
    const volScalarField& T = this->thermo_.T();

    const volScalarField Wm(this->thermo_.W());

    forAll(Y, i)
    {
        const dimensionedScalar Wi
        (
            "W",
            dimMass/dimMoles,
            this->thermo_.composition().Wi(i)
        );
    }

    forAll(Dm_, i)
    {
            PtrList<volScalarField> sumXbyD(Y.size());
            sumXbyD.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "sumXbyD",
                        T.mesh().time().timeName(),
                        T.mesh()
                    ),
                    T.mesh(),
                    dimensionedScalar("sumXbyD", dimensionSet(0,-2,1,0,0,0,0), Foam::SMALL)
                )
            );

            forAll(Y, j)
            {
                if (j != i)
                {
                    sumXbyD[i] += X_[j] / Dij(i,j);
                }
            }

            Dm_[i] = (1-X_[i]) / (sumXbyD[i] + dimensionedScalar("SMALL", dimensionSet(0,-2,1,0,0), Foam::SMALL));
            DmK_[i] = Dm_[i];
    }

    //- Overwrite the values for the porous zones
    forAll(this->pZones_, zoneI)
    {
        const word zoneName = this->pZones_[zoneI].name();
        label znId = this->mesh_.cellZones().findZoneID(zoneName);
        const labelList& cells = this->mesh_.cellZones()[znId];

        const dimensionedScalar RR("RR", dimensionSet(1,2,-2,-1,-1,0,0), 8314.51);
        const scalar PI = 3.14159;

        forAll(cells, cellI)
        {
            label porousID = cells[cellI];

            forAll(DK_, i)
            {
                DK_[i][porousID] = 48.5*dp_[zoneI] / 3* sqrt(8*RR.value()*T[porousID])/(PI*composition.Wi(i));
            }
        }


        forAll(DmK_, i)
        {

            forAll(cells, cellI)
            {
                label porousID = cells[cellI];

                // To prevent dividing by zero due to the fact that Xi of a specie is zero and so Dmi would be, add a SMALL number
                // DmK_[i][porousID] = (eps_[zoneI]/tau_[zoneI]) * (1/( 1/(Dm_[i][porousID]+ Foam::SMALL) + 1/DK_[i][porousID]));    // Still a diffusion coefficient
                DmK_[i][porousID] = (1/( 1/(Dm_[i][porousID]+ Foam::SMALL) + 1/DK_[i][porousID]));    // Still a diffusion coefficient
                DmK_[i][porousID] *= pow(this->phase_[porousID],1.5);
            }
        }
    }
}


template<class Thermo>
void Foam::multiSpeciesTransportModels::FickianBosanquet<Thermo>::updateBinaryDiffusionCoeff()
{
    const PtrList<volScalarField>& Y = this->thermo_.composition().Y();

    for(label i=0; i < Y.size(); i++)
    {
        for(label j=i; j < Y.size(); j++)
        {
            label k = Y.size()*i+j-0.5*i*(i+1);

            Dij_[k] = DijModels_[k].D();
        }
    }
}

template<class Thermo>
const Foam::volScalarField& Foam::multiSpeciesTransportModels::FickianBosanquet<Thermo>::Dij(label i, label j) const
{
    const label species = this->thermo_.composition().Y().size();

    label iStar = min(i,j);
    label jStar = max(i,j);
    label k = species*iStar+jStar-0.5*iStar*(iStar+1);
    return Dij_[k];
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
