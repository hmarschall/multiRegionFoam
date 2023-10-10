/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Fork:     foam-extend
    \\  /    A nd           | Version:  4.1                                 
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

SourceFiles
    genericRegionCoupledJumpFvPatchField.C

Author
    Constantin Habes (constantin.habes@tu-darmstadt.de, Affiliation A)
    Holger Marschall (holger.marschall@tu-darmstadt.de, Affiliation A)

Contact
    Holger Marschall (holger.marschall@tu-darmstadt.de)
    main developer and principal investigator
    TU Darmstadt

Affiliations
    Affiliation A)
    Computational Multiphase Flow
    Department of Mathematics
    Technical University of Darmstadt, Germany

Acknowledgement
    Funded by
    Hessian Ministry of Higher Education, Research, Science and the Arts,
    National High Performance Computing Center for Computational Engineering
    Science (NHR4CES)

Description
    Generic region coupled jump boundary condition for size Psi. Is intended to
    serve as a base class for coupled dirichlet boundary conditions. Includes a
    zero "jump" condition across the interface that can be overridden by child
    classes.

    This file is part of the multiRegionFoam library.

    You may refer to this software as :
    https://doi.org/10.48550/arXiv.2306.01924

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#ifndef genericRegionCoupledJumpFvPatchField_H
#define genericRegionCoupledJumpFvPatchField_H

#include "accelerationModel.H"
#include "interfaceToInterfaceCoupleManager.H"
#include "regionInterfaceType.H"
#include "fixedValueFvPatchFields.H"
#include "gradScheme.H"

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

        //- Neighbour names
        word neighbourRegionName_;

        word neighbourPatchName_;

        word neighbourFieldName_;

protected:

    // Protected data

        //- Diffusivity coefficient
        word kName_;

        //- Jump coefficient
        word KName_;

        //- Acceleration model for BC update
        autoPtr<accelerationModel<Type>> accModel_;

        //- NonOrthogonal correction
        Switch nonOrthCorr_;

        //- Second order accuracy correction
        Switch secondOrder_;

    // Protected member functions

        //- Switch to use directly computed value instead of value jump
        virtual bool useDirectValue() const
        {
            return false;
        }

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
            notImplemented
            (
                "genericRegionCoupledJumpFvPatchField.C\n"
                "tmp<Field<Type>> valueJump()\n"
                "not implemented"
            );

            return (*this * 0);
        };

        //- If direct flux should be computed instead of flux jump
        virtual tmp<Field<Type>> value() const
        {
            notImplemented
            (
                "genericRegionCoupledFluxFvPatchField.C\n"
                "tmp<Field<Type>> fluxJump()\n"
                "not implemented"
            );

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
        scalar normResidual() const;

        //- Return the normalized coupled patch residual
        //- normalised similar to linear system solver residuals
        scalar ofNormResidual() const;

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
