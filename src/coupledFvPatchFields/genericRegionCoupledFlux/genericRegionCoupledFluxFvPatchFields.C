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

#include "objectRegistry.H"
#include "volFields.H"
#include "genericRegionCoupledFluxFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeTemplatePatchTypeField
  (
      fvPatchScalarField,
      genericRegionCoupledFluxFvPatchScalarField
  );

makeTemplatePatchTypeField
  (
      fvPatchVectorField,
      genericRegionCoupledFluxFvPatchVectorField
  );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

template<>
void Foam::genericRegionCoupledFluxFvPatchField<Foam::scalar>::evaluate
(
  const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const volScalarField& vsf =
        this->db().objectRegistry::lookupObject<volScalarField>
        (
            this->dimensionedInternalField().name()
        );

    autoPtr<fvPatchField<vector> > gradpsiPtr(nullptr);

    if
    (
        this->db().objectRegistry::foundObject<volVectorField>
        ("grad(" + this->dimensionedInternalField().name() + ")")
    )
    {
        const fvPatchField<vector>& gradpsi_ =
            patch().lookupPatchField<volVectorField, vector>
            (
                "grad(" + this->dimensionedInternalField().name() + ")"
            );

        gradpsiPtr.reset
        (
            new fvPatchField<vector>(gradpsi_)
        );
    }
    else
    {
        gradpsiPtr.reset
        (
            new fvPatchField<vector>
            (
                fvc::grad(vsf)()
               .boundaryField()[patch().index()]
            )
        );
    }

    const fvPatchField<vector>& gradpsi = gradpsiPtr();

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = delta - n*(n&delta);

    scalarField dpsiP(this->patch().size(), 0);

    if (nonOrthCorr_)
    {
        //- TODO: Add tamplatable nonOrth correction
        dpsiP = (k&gradpsi.patchInternalField());
    }

    if (secondOrder_)
    {
        //- TODO: add tamplatable second order gradient
        scalarField nGradpsiP = (n&gradpsi.patchInternalField());

        Field<scalar>::operator=
        (
            this->patchInternalField() + dpsiP
          + 0.5*(gradient() + nGradpsiP)/this->patch().deltaCoeffs()
        );
    }
    else
    {
            Field<scalar>::operator=
            (
                this->patchInternalField() + dpsiP
              + gradient()/this->patch().deltaCoeffs()
            );
    }

    fvPatchField<scalar>::evaluate();
}

template<>
void Foam::genericRegionCoupledFluxFvPatchField<Foam::vector>::evaluate
(
  const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    fixedGradientFvPatchVectorField::evaluate();

    const volVectorField& vvf =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    autoPtr<fvPatchField<tensor> > gradpsiPtr(nullptr);

    if
    (
        this->db().objectRegistry::foundObject<volTensorField>
        ("grad(" + this->dimensionedInternalField().name() + ")")
    )
    {
        const fvPatchField<tensor>& gradpsi_ =
            patch().lookupPatchField<volTensorField, tensor>
            (
                "grad(" + this->dimensionedInternalField().name() + ")"
            );

        gradpsiPtr.reset
        (
            new fvPatchField<tensor>(gradpsi_)
        );
    }
    else
    {
        gradpsiPtr.reset
        (
            new fvPatchField<tensor>
            (
                fvc::grad(vvf)()
               .boundaryField()[patch().index()]
            )
        );
    }

    const fvPatchField<tensor>& gradpsi = gradpsiPtr();

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = delta - n*(n&delta);

    vectorField dpsiP(this->patch().size(), vector::zero);

    if (nonOrthCorr_)
    {
        //- TODO: Add tamplatable nonOrth correction
        dpsiP = (k&gradpsi.patchInternalField());
    }

    if (secondOrder_)
    {
        //- TODO: add tamplatable second order gradient
        vectorField nGradpsiP = (n&gradpsi.patchInternalField());

        Field<vector>::operator=
        (
            this->patchInternalField() + dpsiP
          + 0.5*(gradient() + nGradpsiP)/this->patch().deltaCoeffs()
        );
    }
    else
    {
            Field<vector>::operator=
            (
                this->patchInternalField() + dpsiP
              + gradient()/this->patch().deltaCoeffs()
            );
    }

    fvPatchField<vector>::evaluate();
}

// ************************************************************************* //
