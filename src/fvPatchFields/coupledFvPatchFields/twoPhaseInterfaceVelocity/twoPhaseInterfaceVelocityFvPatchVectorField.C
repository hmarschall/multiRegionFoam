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

#include "twoPhaseInterfaceVelocityFvPatchVectorField.H"
#include "symmTransformField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

#include "freeSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twoPhaseInterfaceVelocityFvPatchVectorField::
twoPhaseInterfaceVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    ggiFvPatchField<vector>(p, iF),
    master_(false),
    forceJump_(p.size(), vector::zero)
{}


twoPhaseInterfaceVelocityFvPatchVectorField::
twoPhaseInterfaceVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    ggiFvPatchField<vector>(p, iF, dict),
    master_(dict.lookup("master")),
    forceJump_(p.size(), vector::zero)
{
    if (dict.found("forceJump"))
    {
        forceJump_ = vectorField("forceJump", dict, p.size());
    }

//     Info << this->patch().name() << " " << master_ << endl;
}


twoPhaseInterfaceVelocityFvPatchVectorField::
twoPhaseInterfaceVelocityFvPatchVectorField
(
    const twoPhaseInterfaceVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    ggiFvPatchField<vector>(ptf, p, iF, mapper),
    master_(ptf.master_),
    forceJump_(ptf.forceJump_, mapper)
{}


twoPhaseInterfaceVelocityFvPatchVectorField::
twoPhaseInterfaceVelocityFvPatchVectorField
(
    const twoPhaseInterfaceVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    ggiFvPatchField<vector>(ptf, iF),
    master_(ptf.master_),
    forceJump_(ptf.forceJump_)    
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// void twoPhaseInterfaceVelocityFvPatchVectorField::updateCoeffs()
// {
//     if (updated())
//     {
//         return;
//     }
// }

void twoPhaseInterfaceVelocityFvPatchVectorField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    freeSurface& fs =
        const_cast<freeSurface&>
        (
            mesh.lookupObject<freeSurface>("freeSurfaceProperties")
        );

    if (master_)
    {
        scalarField delta = 1.0/this->patch().deltaCoeffs();
        scalarField deltaN = delta*this->patch().weights();
        scalarField deltaP = delta - deltaN;

        vectorField n = this->patch().nf();

        vectorField UP = this->patchInternalField();
        vectorField UN = this->patchNeighbourField();


        // Tangential component
        vectorField UtP = UP - n*(n & UP);
        vectorField UtN = UN - n*(n & UN);

        vectorField Utf = 
            fs.muFluidA().value()*UtP/deltaP
          + fs.muFluidB().value()*UtN/deltaN;

        Utf += (fs.muFluidB().value() - fs.muFluidA().value())
           *(fac::grad(fs.Us()) & fs.aMesh().faceAreaNormals())()
           .internalField();

        // Set tangential surface tension force - to do
        vectorField tangentialSurfaceTensionForce(n.size(), vector::zero);
        Utf += tangentialSurfaceTensionForce;

        Utf /= fs.muFluidA().value()/deltaP
          + fs.muFluidB().value()/deltaN + VSMALL;


        // Normal component
        vectorField UnP = n*(n & UP);
        vectorField UnN = n*(n & UN);

        vectorField Unf = 
            2*fs.muFluidA().value()*UnP/deltaP
          + 2*fs.muFluidB().value()*UnN/deltaN;

        Unf /= 2*fs.muFluidA().value()/deltaP
          + 2*fs.muFluidB().value()/deltaN + VSMALL;


//         UtP -= n*(n & UtP);
//         UtN -= n*(n & UtN);

//         vectorField Utf =
//             fs.muFluidA().value()*UtP/deltaP
//           + fs.muFluidB().value()*UtN/deltaN;

//         vectorField Unf =
//             n*fs.phi().boundaryField()[fs.aPatchID()]
//            /this->patch().magSf();

//         fs.Us().internalField() +=
//             Unf - n*(n & fs.Us().internalField());
//         fs.correctUsBoundaryConditions();
        
//         Utf -= (fs.muFluidA().value() - fs.muFluidB().value())
//            *(fac::grad(fs.Us()) & fs.aMesh().faceAreaNormals())()
//            .internalField();

//         // Set tangential surface tension force - to do
//         vectorField tangentialSurfaceTensionForce(n.size(), vector::zero);

//         Utf += tangentialSurfaceTensionForce;

//         Utf /= fs.muFluidA().value()/deltaP
//           + fs.muFluidB().value()/deltaN + VSMALL;

        vectorField Uf = Utf + Unf;

        fs.Us().internalField() = Uf;
        fs.correctUsBoundaryConditions();

        fs.updateNGradUn();

        Field<vector>::operator=(Uf);
    }
    else
    {
        const ggiFvPatch& ggiPatch = 
            refCast<const ggiFvPatch>(this->patch());

        vectorField Uf = 
            ggiPatch.interpolate
            (
                fs.U().boundaryField()[ggiPatch.shadowIndex()]
            );

        Field<vector>::operator=(Uf);
    }

//     Info << "twoPhaseInterfaceVelocityFvPatchVectorField::initEvaluate: " 
//         << this->patch().name() << endl;
}


// template<class Type>
// void twoPhaseInterfaceVelocityFvPatchVectorField<Type>::evaluate
// (
//     const Pstream::commsTypes
// )
// {
//     fvPatchField<Type>::evaluate();
// }


const Switch& twoPhaseInterfaceVelocityFvPatchVectorField::master() const
{
//     const fvMesh& mesh = this->patch().boundaryMesh().mesh();

//     const freeSurface& fs = 
//         mesh.lookupObject<freeSurface>("freeSurface");

//     if (fs.aPatchID() != this->patch().index())
//     {
//     }

    return master_;
}


void twoPhaseInterfaceVelocityFvPatchVectorField::manipulateMatrix
(
    fvMatrix<vector>& eqn
)
{
    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    const freeSurface& fs =
        mesh.lookupObject<freeSurface>("freeSurfaceProperties");

    vectorField& source = eqn.source();

    vectorField n = this->patch().nf();

    scalarField delta = 1.0/this->patch().deltaCoeffs();
    scalarField deltaN = delta*this->patch().weights();
    scalarField deltaP = delta - deltaN;

    const scalarField& magSf = this->patch().magSf();

    if (master_)
    {
        vectorField tangentialSurfaceTensionForce(n.size(), vector::zero);

//         vectorField newForceJump =
        forceJump_ =
            (fs.muFluidB().value() - fs.muFluidA().value())
           *(fac::grad(fs.Us()) & fs.aMesh().faceAreaNormals())()
           .internalField()
          + tangentialSurfaceTensionForce;

//         forceJump_ += 1.0*(newForceJump - forceJump_); 

        const unallocLabelList& faceCells = this->patch().faceCells();

        forAll(faceCells, faceI)
        {
            source[faceCells[faceI]] += 
                forceJump_[faceI]*magSf[faceI]
               *(fs.muFluidA().value()/deltaP[faceI])
               /(
                    fs.muFluidA().value()/deltaP[faceI]
                  + fs.muFluidB().value()/deltaN[faceI]
                );
        }
    }
    else
    {
        const ggiFvPatch& ggiPatch = 
            refCast<const ggiFvPatch>(this->patch());

        const twoPhaseInterfaceVelocityFvPatchVectorField& shadowPatchU = 
            refCast<const twoPhaseInterfaceVelocityFvPatchVectorField>
            (
                fs.U().boundaryField()[ggiPatch.shadowIndex()]
            );

        forceJump_ = ggiPatch.interpolate(shadowPatchU.forceJump());

        const unallocLabelList& faceCells = this->patch().faceCells();

        forAll(faceCells, faceI)
        {
            source[faceCells[faceI]] += 
                forceJump_[faceI]*magSf[faceI]
               *(fs.muFluidB().value()/deltaP[faceI])
               /(
                    fs.muFluidA().value()/deltaN[faceI]
                  + fs.muFluidB().value()/deltaP[faceI]
                );
        }
    }

    fvPatchField<vector>::manipulateMatrix(eqn);
}


void twoPhaseInterfaceVelocityFvPatchVectorField::patchInterpolate
(
    GeometricField<vector, fvsPatchField, surfaceMesh>& fField,
    const scalarField& pL
) const
{
    fField.boundaryField()[this->patch().index()] = *this;

//     fField.boundaryField()[patchI] =
//         pL*this->patchInternalField()
//       + (1 - pL)*this->patchNeighbourField();
}


void twoPhaseInterfaceVelocityFvPatchVectorField::patchInterpolate
(
    GeometricField<vector, fvsPatchField, surfaceMesh>& fField,
    const scalarField& pL,
    const scalarField& pY
) const
{
    fField.boundaryField()[this->patch().index()] = *this;

//     fField.boundaryField()[patchI] =
//         pL*this->patchInternalField()
//       + pY*this->patchNeighbourField();
}


void twoPhaseInterfaceVelocityFvPatchVectorField::
write(Ostream& os) const
{
    ggiFvPatchField<vector>::write(os);

    os.writeKeyword("master")
        << master_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    twoPhaseInterfaceVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
