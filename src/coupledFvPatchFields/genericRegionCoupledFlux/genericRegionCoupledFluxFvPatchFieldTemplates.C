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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    genericRegionCoupledFluxFvPatchFieldFvPatchFieldTemplates

Description
    Template specialisations

SourceFiles
    genericRegionCoupledFluxFvPatchFieldFvPatchFieldTemplates.C


\*---------------------------------------------------------------------------*/

#include "genericRegionCoupledFluxFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void genericRegionCoupledFluxFvPatchField<scalar>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const fvPatchField<vector>& gradpsi =
    patch().lookupPatchField<volVectorField, vector>
    (
        "grad(" + psiName_ + ")"
    );

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

    updatePhi();

    fvPatchField<scalar>::evaluate();
}

template<>
void genericRegionCoupledFluxFvPatchField<vector>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    fixedGradientFvPatchVectorField::evaluate();

    const fvPatchField<tensor>& gradpsi =
    patch().lookupPatchField<volTensorField, tensor>
    (
        "grad(" + psiName_ + ")"
    );

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

    updatePhi();

    fvPatchField<vector>::evaluate();
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
