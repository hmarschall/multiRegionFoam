/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "heRhoThermo.H"
#include "fvMesh.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// template<class BasicThermo, class MixtureType>
// void Foam::heRhoThermo<BasicThermo, MixtureType>::calculate()
// {
//     const scalarField& hCells = this->he_.internalField();
//     const scalarField& pCells = this->p_.internalField();

//     scalarField& TCells = this->T_.internalField();
//     scalarField& psiCells = this->psi_.internalField();
//     scalarField& rhoCells = this->rho_.internalField();
//     scalarField& muCells = this->mu_.internalField();
//     scalarField& alphaCells = this->alpha_.internalField();

//     forAll(TCells, celli)
//     {
//         const typename MixtureType::thermoType& mixture =
//             this->cellMixture(celli);

//         TCells[celli] = mixture.THE(hCells[celli],pCells[celli], TCells[celli]);  // THE which is actually returning the mixture of specie values THs(he, p, T0); -> do not use sensibleEnthalpy e.g. here 
//         psiCells[celli] = mixture.psi(pCells[celli], TCells[celli]);
//         rhoCells[celli] = mixture.rho(pCells[celli], TCells[celli]);

//         muCells[celli] = mixture.mu(pCells[celli], TCells[celli]);
//         alphaCells[celli] = mixture.alphah(pCells[celli], TCells[celli]);
//     }

//     forAll(this->T_.boundaryField(), patchi)
//     {
//         fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
//         fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
//         fvPatchScalarField& ppsi = this->psi_.boundaryField()[patchi];
//         fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];

//         fvPatchScalarField& ph = this->he_.boundaryField()[patchi];

//         fvPatchScalarField& pmu_ = this->mu_.boundaryField()[patchi];
//         fvPatchScalarField& palpha_ = this->alpha_.boundaryField()[patchi];

//         if (pT.fixesValue())
//         {
//             forAll(pT, facei)
//             {
//                 const typename MixtureType::thermoType& mixture =
//                     this->patchFaceMixture(patchi, facei);

//                 ph[facei] = mixture.H(pT[facei]);

//                 ppsi[facei] = mixture.psi(pp[facei], pT[facei]);
//                 prho[facei] = mixture.rho(pp[facei], pT[facei]);
//                 pmu_[facei] = mixture.mu(pT[facei]);
//                 palpha_[facei] = mixture.alphah(pT[facei]);
//             }
//         }
//         else
//         {
//             forAll(pT, facei)
//             {
//                 const typename MixtureType::thermoType& mixture =
//                     this->patchFaceMixture(patchi, facei);

//                 pT[facei] = mixture.THE(ph[facei], pT[facei]);

//                 ppsi[facei] = mixture.psi(pp[facei], pT[facei]);
//                 prho[facei] = mixture.rho(pp[facei], pT[facei]);
//                 pmu_[facei] = mixture.mu(pT[facei]);
//                 palpha_[facei] = mixture.alphah(pT[facei]);
//             }
//         }
//     }
// }


template<class BasicThermo, class MixtureType>
void Foam::heRhoThermo<BasicThermo, MixtureType>::calculate
(
    const volScalarField& p,
    volScalarField& T,
    volScalarField& he,
    volScalarField& psi,
    volScalarField& rho,
    volScalarField& mu,
    volScalarField& alpha,
    const bool doOldTimes
)
{
    // Note: update oldTimes before current time so that if T.oldTime() is
    // created from T, it starts from the unconverted T
    if (doOldTimes && (p.nOldTimes() || T.nOldTimes()))
    {
        calculate
        (
            p.oldTime(),
            T.oldTime(),
            he.oldTime(),
            psi.oldTime(),
            rho.oldTime(),
            mu.oldTime(),
            alpha.oldTime(),
            true
        );
    }

    const scalarField& hCells = he.internalField();
    const scalarField& pCells = p.internalField();

    scalarField& TCells = T.internalField();
    scalarField& psiCells = psi.internalField();
    scalarField& rhoCells = rho.internalField();
    scalarField& muCells = mu.internalField();
    scalarField& alphaCells = alpha.internalField();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        if (this->updateT())
        {
            TCells[celli] = mixture_.THE
            (
                hCells[celli],
                pCells[celli],
                TCells[celli]
            );
        }

        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& ppsi = this->psi_.boundaryField()[patchi];
        fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        fvPatchScalarField& phe = this->he_.boundaryField()[patchi];
        fvPatchScalarField& pmu = this->mu_.boundaryField()[patchi];
        fvPatchScalarField& palpha = this->alpha_.boundaryField()[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                if (this->updateT())
                {
                    pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);
                }

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
Foam::heRhoThermo<BasicThermo, MixtureType>::heRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName           // obj is the mesh (see hREactionThermo.C)
)
:
    heThermo<BasicThermo, MixtureType>(mesh, phaseName)
{
    // // Info << "Print constructor heRhoThermo 1 " << endl;
    // // scalarField& hCells = this->he_.internalField();
    // // const scalarField& TCells = this->T_.internalField();

    // // Info << "Print constructor heRhoThermo 2 " << endl;

    // // // forAll(hCells, celli)
    // // // {
    // // //     hCells[celli] = this->cellMixture(celli).H(TCells[celli]);
    // // // }

    // // Info << "Print constructor heRhoThermo 3 " << endl;

    // // // forAll(this->h_.boundaryField(), patchi)
    // // // {
    // // //     this->h_.boundaryField()[patchi] == h(this->T_.boundaryField()[patchi], patchi);
    // // // }
    // // Info << "Print constructor heRhoThermo 4 " << endl;

    // // // this->hBoundaryCorrection(this->h_);

    // // // calculate();
    // // Info << "Print constructor heRhoThermo 5 " << endl;

    calculate
    (
        this->p_,
        this->T_,
        this->he_,
        this->psi_,
        this->rho_,
        this->mu_,
        this->alpha_,
        true                    // Create old time fields
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
Foam::heRhoThermo<BasicThermo, MixtureType>::~heRhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
void Foam::heRhoThermo<BasicThermo, MixtureType>::correct()
{
    if (debug)
    {
        Info<< "entering heRhoThermo<BasicThermo, MixtureType>::correct()" << endl;
    }

    // calculate();
    calculate
    (
        this->p_,
        this->T_,
        this->he_,
        this->psi_,
        this->rho_,
        this->mu_,
        this->alpha_,
        false           // No need to update old times
    );

    if (debug)
    {
        Info<< "exiting heRhoThermo<BasicThermo, MixtureType>::correct()" << endl;
    }
}

// ************************************************************************* //
