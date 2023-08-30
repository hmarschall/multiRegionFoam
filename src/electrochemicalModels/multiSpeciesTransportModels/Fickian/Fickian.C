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

#include "Fickian.H"
#include "fvmLaplacian.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::multiSpeciesTransportModels::Fickian<Thermo>::Fickian
(
    const phaseModel& phase
)
:
    MultiSpeciesTransportModel<Thermo>(phase),
    X_(phase.X()),
    Dm_(this->thermo_.composition().species().size()),
    DijModels_(0.5*this->thermo_.composition().species().size()*(this->thermo_.composition().species().size()+1)),
    Dij_(DijModels_.size())
{
    const PtrList<volScalarField>& Y = phase.Y();

    forAll(Y, i)
    {
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
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::multiSpeciesTransportModels::Fickian<Thermo>::~Fickian()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::tmp<Foam::volScalarField> Foam::multiSpeciesTransportModels::Fickian<Thermo>::DEff
(
    const volScalarField& Yi
) const
{
    const basicSpecieMixture& composition = this->thermo_.composition();

    return this->thermo_.rho()*Dm_[composition.index(Yi)];
}


template<class Thermo>
Foam::tmp<Foam::scalarField> Foam::multiSpeciesTransportModels::Fickian<Thermo>::DEff
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
Foam::tmp<Foam::surfaceScalarField> Foam::multiSpeciesTransportModels::Fickian<Thermo>::q() const
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
Foam::tmp<Foam::fvScalarMatrix> Foam::multiSpeciesTransportModels::Fickian<Thermo>::divq
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
Foam::tmp<Foam::surfaceScalarField> Foam::multiSpeciesTransportModels::Fickian<Thermo>::j
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

    tj = -fvc::interpolate(this->phase_*DEff(Yi))*fvc::snGrad(Yi);      // This would call Fickian:.DEff (if part of template) 

    return tj;
}


template<class Thermo>
Foam::tmp<Foam::fvScalarMatrix> Foam::multiSpeciesTransportModels::Fickian<Thermo>::divj
(
    volScalarField& Yi
) const
{
    return -fvm::laplacian(this->phase_*DEff(Yi), Yi);
}


template<class Thermo>
void Foam::multiSpeciesTransportModels::Fickian<Thermo>::correct()
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

        // X_[i] = Wm*Y[i]/Wi;
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
    }
    // }
}

template<class Thermo>
void Foam::multiSpeciesTransportModels::Fickian<Thermo>::updateBinaryDiffusionCoeff()
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
const Foam::volScalarField& Foam::multiSpeciesTransportModels::Fickian<Thermo>::Dij(label i, label j) const
{
    const label species = this->thermo_.composition().Y().size();

    label iStar = min(i,j);
    label jStar = max(i,j);
    label k = species*iStar+jStar-0.5*iStar*(iStar+1);
    return Dij_[k];
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
