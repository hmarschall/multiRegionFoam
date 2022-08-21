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
    genericRegionCoupledJumpFvPatchFieldFvPatchFieldTemplates

Description
    Template specialisations

SourceFiles
    genericRegionCoupledJumpFvPatchFieldFvPatchFieldTemplates.C


\*---------------------------------------------------------------------------*/

#include "genericRegionCoupledJumpFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
tmp<Field<scalar> > genericRegionCoupledJumpFvPatchField<scalar>::snGrad() const
{

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

        tmp<scalarField> tnGradpsi
        (
            new scalarField(this->patch().size(), 0)
        );

        tnGradpsi() =
            2
           *(
                *this
              - (patchInternalField() + dpsiP)
            )*this->patch().deltaCoeffs()
          - nGradpsiP;

        return tnGradpsi;
    }

    // First order
    tmp<scalarField> tnGradpsi
    (
        new scalarField(this->patch().size(), 0)
    );

    tnGradpsi() =
        (
            *this
          - (patchInternalField() + dpsiP)
        )*this->patch().deltaCoeffs();

    return tnGradpsi;
}

template<>
tmp<Field<vector> > genericRegionCoupledJumpFvPatchField<vector>::snGrad() const
{

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
        dpsiP = (k&gradpsi.patchInternalField());
    }

    if (secondOrder_)
    {
        vectorField nGradpsiP = (n&gradpsi.patchInternalField());

        tmp<vectorField> tnGradpsi
        (
            new vectorField(this->patch().size(), vector::zero)
        );

        tnGradpsi() =
            2
           *(
                *this
              - (patchInternalField() + dpsiP)
            )*this->patch().deltaCoeffs()
          - nGradpsiP;

        tnGradpsi() -= n*(n&tnGradpsi());

        return tnGradpsi;
    }

    // First order
    tmp<vectorField> tnGradpsi
    (
        new vectorField(this->patch().size(), vector::zero)
    );

    tnGradpsi() =
        (
            *this
          - (patchInternalField() + dpsiP)
        )*this->patch().deltaCoeffs();

    tnGradpsi() -= n*(n&tnGradpsi());

    return tnGradpsi;
}

template<>
tmp<Field<scalar> > genericRegionCoupledJumpFvPatchField<scalar>::gradientBoundaryCoeffs() const
{
    const fvPatchField<vector>& gradpsi =
        patch().lookupPatchField<GeometricField<vector, fvPatchField, volMesh>, vector>
        (
            "grad(" + psiName_ + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = delta - n*(n&delta);

    scalarField dpsiP(this->patch().size(), 0);

    if (nonOrthCorr_)
    {
        dpsiP = (k&gradpsi.patchInternalField());
    }

    if (secondOrder_)
    {
        scalarField nGradpsiP = (n&gradpsi.patchInternalField());

        return
            this->patch().deltaCoeffs()
           *(
                2*(*this - dpsiP) 
              - this->patchInternalField()
            )
          - nGradpsiP;
    }

    return this->patch().deltaCoeffs()*(*this - dpsiP);
}

template<>
tmp<Field<vector> > genericRegionCoupledJumpFvPatchField<vector>::gradientBoundaryCoeffs() const
{
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
        dpsiP = (k&gradpsi.patchInternalField());
    }

    if (secondOrder_)
    {
        vectorField nGradpsiP = (n&gradpsi.patchInternalField());

        vectorField nGradpsi =
            2
           *(
                *this
              - (patchInternalField() + dpsiP)
            )*this->patch().deltaCoeffs()
          - nGradpsiP;

        vectorField nGradpsin = n*(n&nGradpsi);

        return
            this->patch().deltaCoeffs()
           *(
                2*(*this - dpsiP) 
              - this->patchInternalField()
            )
          - nGradpsiP
          - nGradpsin;
    }

    // First order
    vectorField nGradpsi =
        (
            *this
          - (this->patchInternalField() + dpsiP)
        )*this->patch().deltaCoeffs();

    vectorField nGradpsin = n*(n&nGradpsi);

    return 
        this->patch().deltaCoeffs()
       *(
            *this - dpsiP
        )
      - nGradpsin;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
