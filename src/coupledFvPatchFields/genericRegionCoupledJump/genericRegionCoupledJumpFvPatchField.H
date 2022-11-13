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

Class
    genericRegionCoupledJumpFvPatchField

Description
    Generic region coupled jump boundary condition for size Psi. Is intended to
    serve as a base class for coupled dirichlet boundary conditions. Includes a 
    zero "jump" condition across the interface that can be overridden by child 
    classes.

    Example usage:
        interfacePatch
        {
            type                    genericRegionCoupledJump;
            k                       "mu";
            relax                   0.5;
            nonOrthCorr             yes;
            secondOrder             yes;
            value                   uniform 0;
        }

SourceFiles
    genericRegionCoupledJumpFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef genericRegionCoupledJumpFvPatchField_H
#define genericRegionCoupledJumpFvPatchField_H

#include "interfaceToInterfaceCoupleManager.H"
#include "regionInterface.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class genericRegionCoupledJumpFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class genericRegionCoupledJumpFvPatchField
:
    public fixedValueFvPatchField<Type>,
    public interfaceToInterfaceCoupleManager
{
    // Private data

protected: 

    //- Diffusion coefficient
    word kName_;

    //- Relaxation factor for BC update
    scalar relax_;

    //- NonOrthogonal correction
    Switch nonOrthCorr_;

    //- Second order accuracy correction
    Switch secondOrder_;

    //- Update phi 
    virtual void updatePhi()
    {
        //- Override in child class if needed
    };

public:

    //- Runtime type information
    TypeName("genericRegionCoupledJump");


    // Constructors

        //- Construct from patch and internal field
        genericRegionCoupledJumpFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        genericRegionCoupledJumpFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given genericRegionCoupledJumpFvPatchField
        //  onto a new patch
        genericRegionCoupledJumpFvPatchField
        (
            const genericRegionCoupledJumpFvPatchField&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new genericRegionCoupledJumpFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        genericRegionCoupledJumpFvPatchField
        (
            const genericRegionCoupledJumpFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new genericRegionCoupledJumpFvPatchField(*this, iF)
            );
        }

    // Member functions

        //- Calculate jump
        virtual tmp<Field<Type>> valueJump() const
        {
                return (*this * 0);
        };

        //- Update the patch field coefficients
        virtual void updateCoeffs();

        //- Return the patch flux
        virtual tmp<Field<Type>> flux() const;

        //- Return patch-normal gradient
        virtual tmp<Field<Type> > snGrad() const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<Type> > gradientBoundaryCoeffs() const;

        //- Return the raw coupled patch residual
        scalarField rawResidual() const;

        //- Return the normalized coupled patch residual
        scalarField normResidual() const;

        //- Return the normalized coupled patch residual
        //- normalised similar to linear system solver residuals
        scalar ofNormResidual() const;

        //- Return the maximum normalized coupled patch residual
        scalar maxNormResidual() const;

        //- Write
        virtual void write(Ostream&) const;
};

//- Specialisations
template<> tmp<Field<scalar> >
genericRegionCoupledJumpFvPatchField<scalar>::snGrad() const;

template<> tmp<Field<vector> >
genericRegionCoupledJumpFvPatchField<vector>::snGrad() const;

template<> tmp<Field<scalar> > 
genericRegionCoupledJumpFvPatchField<scalar>::gradientBoundaryCoeffs() const;

template<> tmp<Field<vector> >
genericRegionCoupledJumpFvPatchField<vector>::gradientBoundaryCoeffs() const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "genericRegionCoupledJumpFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
