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

#include "unityLewisFourier.H"
#include "fvmLaplacian.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::multiSpeciesTransportModels::unityLewisFourier<Thermo>::unityLewisFourier
(
    const phaseModel& phase
)
:
    MultiSpeciesTransportModel<Thermo>(phase),
    X_(phase.X()),
    Dm_(this->thermo_.composition().species().size())
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
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::multiSpeciesTransportModels::unityLewisFourier<Thermo>::~unityLewisFourier()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::tmp<Foam::volScalarField> Foam::multiSpeciesTransportModels::unityLewisFourier<Thermo>::DEff
(
    const volScalarField& Yi
) const
{
    const basicSpecieMixture& composition = this->thermo_.composition();

    Info << "Print Yi.member() " << Yi.member() << endl;
    Info << "Print Dm_[composition.index(Yi)] " << gAverage(Dm_[composition.index(Yi)]) << endl;

    return this->thermo_.kappa()/this->thermo_.Cp();
}


template<class Thermo>
Foam::tmp<Foam::scalarField> Foam::multiSpeciesTransportModels::unityLewisFourier<Thermo>::DEff
(
    const volScalarField& Yi,
    const label patchi
) const
{
    return
         this->thermo_.kappa()().boundaryField()[patchi]
               /this->thermo_.Cp()().boundaryField()[patchi];
}


template<class Thermo>
Foam::tmp<Foam::surfaceScalarField> Foam::multiSpeciesTransportModels::unityLewisFourier<Thermo>::q() const
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
            -fvc::interpolate(this->phase_*this->thermo_.alphahe())
            *fvc::snGrad(this->thermo_.he())
        )
    );

    return tmpq;
}


template<class Thermo>
Foam::tmp<Foam::fvScalarMatrix> Foam::multiSpeciesTransportModels::unityLewisFourier<Thermo>::divq
(
    volScalarField& he
) const
{
   return -fvm::laplacian(this->phase_*this->thermo_.alphahe(), he);
}


template<class Thermo>
Foam::tmp<Foam::surfaceScalarField> Foam::multiSpeciesTransportModels::unityLewisFourier<Thermo>::j
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
Foam::tmp<Foam::fvScalarMatrix> Foam::multiSpeciesTransportModels::unityLewisFourier<Thermo>::divj
(
    volScalarField& Yi
) const
{
    return -fvm::laplacian(this->phase_*DEff(Yi), Yi);
}


template<class Thermo>
void Foam::multiSpeciesTransportModels::unityLewisFourier<Thermo>::correct()
{
    //- Dummy
}

template<class Thermo>
void Foam::multiSpeciesTransportModels::unityLewisFourier<Thermo>::updateBinaryDiffusionCoeff()
{
    //- Dummy
}

template<class Thermo>
const Foam::volScalarField& Foam::multiSpeciesTransportModels::unityLewisFourier<Thermo>::Dij(label i, label j) const
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
