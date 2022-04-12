/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "twoPhaseInterfacePressureFvPatchScalarField.H"
#include "symmTransformField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

#include "freeSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twoPhaseInterfacePressureFvPatchScalarField::
twoPhaseInterfacePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ggiFvPatchField<scalar>(p, iF),
    master_(false),
    jump_(p.size(), 0)
{}


twoPhaseInterfacePressureFvPatchScalarField::
twoPhaseInterfacePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    ggiFvPatchField<scalar>(p, iF, dict),
    master_(dict.lookup("master")),
    jump_(p.size(), 0)
    // jump_("jump", dict, p.size())
{
//     Info << this->patch().name() << " " << master_ << endl;

    if (dict.found("jump"))
    {
        jump_ = scalarField("jump", dict, p.size());
    }

}


twoPhaseInterfacePressureFvPatchScalarField::
twoPhaseInterfacePressureFvPatchScalarField
(
    const twoPhaseInterfacePressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    ggiFvPatchField<scalar>(ptf, p, iF, mapper),
    master_(ptf.master_),
    jump_(ptf.jump_)
{}


twoPhaseInterfacePressureFvPatchScalarField::
twoPhaseInterfacePressureFvPatchScalarField
(
    const twoPhaseInterfacePressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ggiFvPatchField<scalar>(ptf, iF),
    master_(ptf.master_),
    jump_(ptf.jump_)    
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<scalarField> twoPhaseInterfacePressureFvPatchScalarField::curJump() const
{
    tmp<scalarField> tJump(new scalarField(this->patch().size(), 0));

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    const freeSurface& fs = 
        mesh.lookupObject<freeSurface>("freeSurfaceProperties");

    const scalarField& K = fs.aMesh().faceCurvatures().internalField();

    tJump() = 
        fs.cleanInterfaceSurfTension().value()*K
      - 2.0*(fs.muFluidA().value() - fs.muFluidB().value())*fs.nGradUn()
      + (fs.rhoFluidA().value() - fs.rhoFluidB().value())
       *(
            fs.g().value()
          & (
                fs.mesh().C().boundaryField()[fs.aPatchID()]
            )
        );

    return tJump;
}


tmp<Field<scalar> > 
twoPhaseInterfacePressureFvPatchScalarField::patchNeighbourField() const
{
    tmp<Field<scalar> > tpnf
    (
        new Field<scalar>(ggiFvPatchField<scalar>::patchNeighbourField())
    );
    Field<scalar>& pnf = tpnf();

    if (master_)
    {
        for (label faceI = 0; faceI < pnf.size(); faceI++)
        {
            pnf[faceI] -= jump_[faceI];
        }
    }
    else
    {
        for (label faceI = 0; faceI < pnf.size(); faceI++)
        {
            pnf[faceI] += jump_[faceI];
        }        
    }

    return tpnf;
}


void twoPhaseInterfacePressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (master_)
    {
        jump_ = curJump();
    }
    else
    {
        const fvMesh& mesh = this->patch().boundaryMesh().mesh();

        const freeSurface& fs = 
            mesh.lookupObject<freeSurface>("freeSurfaceProperties");

        const ggiFvPatch& ggiPatch = 
            refCast<const ggiFvPatch>(this->patch());

        const twoPhaseInterfacePressureFvPatchScalarField& shadowPatchP = 
            refCast<const twoPhaseInterfacePressureFvPatchScalarField>
            (
                fs.p().boundaryField()[ggiPatch.shadowIndex()]
            );

        jump_ = ggiPatch.interpolate(shadowPatchP.jump());
    }
}


void twoPhaseInterfacePressureFvPatchScalarField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    const freeSurface& fs =
        mesh.lookupObject<freeSurface>("freeSurfaceProperties");

    const ggiFvPatch& ggiPatch = 
        refCast<const ggiFvPatch>(this->patch());

    scalarField pf(this->patch().size(), 0);

    if (master_)
    {
        const unallocLabelList& sfc = ggiPatch.shadow().faceCells();

        // Reconstruct AU
        const volScalarField& AU = 
            mesh.lookupObject<volScalarField>("AU");
        const scalarField& AUI = AU.internalField();

        scalarField sAU(sfc.size(), 0);

        forAll (sAU, i)
        {
            sAU[i] = AUI[sfc[i]];
        }

        scalarField nrAU = 
            1.0/ggiPatch.interpolate(sAU);
        scalarField orAU = 
            1.0/AU.boundaryField()[this->patch().index()].patchInternalField();

        // Reconstruct phi
        scalarField nHoAn = 
            ggiPatch.interpolate
            (
                fs.phi().boundaryField()[ggiPatch.shadow().index()]
               /ggiPatch.shadow().magSf()
            );
        nHoAn *= -1;
        scalarField oHoAn = 
            fs.phi().boundaryField()[this->patch().index()]
           /this->patch().magSf();

        scalarField dHoAn = oHoAn - nHoAn;

        scalarField delta = 1.0/this->patch().deltaCoeffs();
        scalarField nDelta = delta*this->patch().weights();
        scalarField oDelta = delta - nDelta;

        scalarField oP = this->patchInternalField();
        scalarField nP = this->patchNeighbourField();

        pf = orAU*oP/oDelta + nrAU*nP/nDelta + dHoAn;
        pf /= orAU/oDelta + nrAU/nDelta + SMALL;
    }
    else
    {
        pf = 
            ggiPatch.interpolate
            (
                fs.p().boundaryField()[ggiPatch.shadowIndex()]
            )
          + jump_;
    }

    Field<scalar>::operator=(pf);

//     Field<scalar> pf
//     (
//         this->patch().weights()*this->patchInternalField()
//       + (1.0 - this->patch().weights())*this->patchNeighbourField()
//     );

//     Info << this->patch().name() << " dHoAn , min: " << min(dHoAn) 
//         << ", max: " << max(dHoAn) << ", avg: " << average(dHoAn) << endl;

//     Info << this->patch().name() << " pf , min: " << min(pf) 
//         << ", max: " << max(pf) << ", avg: " << average(pf) << endl;
}


// template<class Type>
// void twoPhaseInterfacePressureFvPatchScalarField<Type>::evaluate
// (
//     const Pstream::commsTypes
// )
// {
//     fvPatchField<Type>::evaluate();
// }


void twoPhaseInterfacePressureFvPatchScalarField::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    // Communication is allowed either before or after processor
    // patch comms.  HJ, 11/Jul/2011

    const ggiFvPatch& ggiPatch = refCast<const ggiFvPatch>(this->patch());

    // Get shadow face-cells and assemble shadow field
    const unallocLabelList& sfc = ggiPatch.shadow().faceCells();

    scalarField sField(sfc.size());

    forAll (sField, i)
    {
        sField[i] = psiInternal[sfc[i]];
    }

    scalarField pnf = ggiPatch.interpolate(sField);

//     Info << psiInternal.size() << ", " << this->internalField().size() 
//         << endl;

    if
    (
        reinterpret_cast<const void*>(&psiInternal)
     == reinterpret_cast<const void*>(&this->internalField())
    )
    {
//         Info << "test-1 " << cmpt << endl;
        if (master_)
        {
            for (label faceI = 0; faceI < pnf.size(); faceI++)
            {
                pnf[faceI] -= jump_[faceI];
            }
        }
        else
        {
            for (label faceI = 0; faceI < pnf.size(); faceI++)
            {
                pnf[faceI] += jump_[faceI];
            }
        }
    }
//     else
//     {
//         Info << "test-2 " << cmpt << endl;
//     }

    // Multiply the field by coefficients and add into the result
    const unallocLabelList& fc = ggiPatch.faceCells();

    if (switchToLhs)
    {
        forAll(fc, elemI)
        {
            result[fc[elemI]] += coeffs[elemI]*pnf[elemI];
        }
    }
    else
    {
        forAll(fc, elemI)
        {
            result[fc[elemI]] -= coeffs[elemI]*pnf[elemI];
        }
    }
}

const Switch& twoPhaseInterfacePressureFvPatchScalarField::master() const
{
//     const fvMesh& mesh = this->patch().boundaryMesh().mesh();

//     const freeSurface& fs = 
//         mesh.lookupObject<freeSurface>("freeSurface");

//     if (fs.aPatchID() != this->patch().index())
//     {
//     }

    return master_;
}


void twoPhaseInterfacePressureFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& eqn
)
{
    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    const freeSurface& fs =
        mesh.lookupObject<freeSurface>("freeSurfaceProperties");

    const ggiFvPatch& ggiPatch = 
        refCast<const ggiFvPatch>(this->patch());

    const unallocLabelList& sfc = ggiPatch.shadow().faceCells();

    // Reconstruct AU
    const volScalarField& AU = 
        mesh.lookupObject<volScalarField>("AU");
    const scalarField& AUI = AU.internalField();

    scalarField sAU(sfc.size(), 0);

    forAll (sAU, i)
    {
        sAU[i] = AUI[sfc[i]];
    }

    scalarField nrAU = 
        1.0/ggiPatch.interpolate(sAU);
    scalarField orAU = 
        1.0/AU.boundaryField()[this->patch().index()].patchInternalField();

    // Reconstruct phi
    scalarField nHoAn = 
        ggiPatch.interpolate
        (
            fs.phi().boundaryField()[ggiPatch.shadow().index()]
           /ggiPatch.shadow().magSf()
        );
    nHoAn *= -1;
    scalarField oHoAn = 
        fs.phi().boundaryField()[this->patch().index()]
       /this->patch().magSf();

    scalarField dHoAn(ggiPatch.size(), 0);
    if (master_)
    {
        dHoAn = oHoAn - nHoAn;
    }
    else
    {
        dHoAn = -(nHoAn - oHoAn);
    }

//     dHoAn *= 0;
//     vectorField n = this->patch().nf();

    scalarField delta = 1.0/this->patch().deltaCoeffs();
    scalarField nDelta = delta*this->patch().weights();
    scalarField oDelta = delta - nDelta;

    scalarField pSource = 
        dHoAn
       *(orAU/oDelta)
       /(
            orAU/oDelta
          + nrAU/nDelta
        );

    const unallocLabelList& faceCells = this->patch().faceCells();

    const scalarField& magSf = this->patch().magSf();

    scalarField& source = eqn.source();

    forAll(faceCells, faceI)
    {
        source[faceCells[faceI]] -= pSource[faceI]*magSf[faceI];
    }

//     if (master_)
//     {
//         forAll(faceCells, faceI)
//         {
//             source[faceCells[faceI]] -= pSource[faceI]*magSf[faceI];
//         }
//     }
//     else
//     {
//         forAll(faceCells, faceI)
//         {
//             source[faceCells[faceI]] += pSource[faceI]*magSf[faceI];
//         }
//     }

    fvPatchField<scalar>::manipulateMatrix(eqn);
}


void twoPhaseInterfacePressureFvPatchScalarField::patchInterpolate
(
    GeometricField<scalar, fvsPatchField, surfaceMesh>& fField,
    const scalarField& pL
) const
{
    
    fField.boundaryField()[this->patch().index()] = *this;

//     Info << "_" << this->patch().name() << ", min: " << min(*this) 
//         << ", max: " << max(*this) 
//         << ", avg: " << average(*this) << endl;

//     fField.boundaryField()[patchI] =
//         pL*this->patchInternalField()
//       + (1 - pL)*this->patchNeighbourField();
}


void twoPhaseInterfacePressureFvPatchScalarField::patchInterpolate
(
    GeometricField<scalar, fvsPatchField, surfaceMesh>& fField,
    const scalarField& pL,
    const scalarField& pY
) const
{
    fField.boundaryField()[this->patch().index()] = *this;

//     fField.boundaryField()[patchI] =
//         pL*this->patchInternalField()
//       + pY*this->patchNeighbourField();
}


void twoPhaseInterfacePressureFvPatchScalarField::patchFlux
(
    GeometricField<scalar, fvsPatchField, surfaceMesh>& flux,
    const fvMatrix<scalar>& matrix
) const
{
    const label patchI = this->patch().index();

//         // Coupled patch
//         flux.boundaryField()[patchI] =
//             cmptMultiply
//             (
//                 matrix.internalCoeffs()[patchI],
//                 this->patchInternalField()
//             )
//           - cmptMultiply
//             (
//                 matrix.boundaryCoeffs()[patchI],
//                 this->patchNeighbourField()
//             );

    const fvMesh& mesh = 
        this->patch().boundaryMesh().mesh();

    const volScalarField& AU = 
        mesh.lookupObject<volScalarField>("AU");

    scalarField orAU = 
        1.0/AU.boundaryField()[patchI].patchInternalField();

    scalarField delta = 1.0/this->patch().deltaCoeffs();
    scalarField nDelta = delta*this->patch().weights();
    scalarField oDelta = delta - nDelta;

    flux.boundaryField()[patchI] = 
        this->patch().magSf()
       *orAU*(*this - this->patchInternalField())/oDelta;

//     Info << this->patch().name() << " " 
//         << sum(mag(flux.boundaryField()[patchI])) << endl;
}

void twoPhaseInterfacePressureFvPatchScalarField::write(Ostream& os) const
{
    ggiFvPatchField<scalar>::write(os);

    os.writeKeyword("master")
        << master_ << token::END_STATEMENT << nl;

    jump_.writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    twoPhaseInterfacePressureFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
