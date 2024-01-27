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
\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "primitiveFieldsFwd.H"
namespace Foam
{

// Forward declaration of generic coupled BCs
template<class Type> class genericRegionCoupledJumpFvPatchField;
template<class Type> class genericRegionCoupledFluxFvPatchField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > regionInterfaceType::interpolateFacesFromA
(
    const Field<Type>& fromField
) const
{
    notImplemented
    (
        "regionInterfaceTypeTemplates.C\n"
        "tmp<Field<Type> > regionInterfaceType::interpolateFacesFromA\n"
        "(\n"
        "const Field<Type>& fromField\n"
        ") const\n"
        "not implemented"
    );

    tmp<Field<Type> > ttoField
    (
        new Field<Type>(globalPatchB().globalPatch().size())
    );

    return ttoField;
}

template<class Type>
tmp<Field<Type> > regionInterfaceType::interpolateFacesFromB
(
    const Field<Type>& fromField
) const
{
    notImplemented
    (
        "regionInterfaceTypeTemplates.C\n"
        "tmp<Field<Type> > regionInterfaceType::interpolateFacesFromB\n"
        "(\n"
        "const Field<Type>& fromField\n"
        ") const\n"
        "not implemented"
    );

    tmp<Field<Type> > ttoField
    (
        new Field<Type>(globalPatchA().globalPatch().size())
    );

    return ttoField;
}

template<class Type>
tmp<Field<Type> > regionInterfaceType::interpolatePointsFromA
(
    const Field<Type>& fromField
) const
{
    notImplemented
    (
        "regionInterfaceTypeTemplates.C\n"
        "tmp<Field<Type> > regionInterfaceType::interpolatePointsFromA\n"
        "(\n"
        "const Field<Type>& fromField\n"
        ") const\n"
        "not implemented"
    );

    tmp<Field<Type> > ttoField
    (
        new Field<Type>(globalPatchB().globalPatch().size())
    );

    return ttoField;
}

template<class Type>
tmp<Field<Type> > regionInterfaceType::interpolatePointsFromB
(
    const Field<Type>& fromField
) const
{
    notImplemented
    (
        "regionInterfaceTypeTemplates.C\n"
        "tmp<Field<Type> > regionInterfaceType::interpolatePointsFromB\n"
        "(\n"
        "const Field<Type>& fromField\n"
        ") const\n"
        "not implemented"
    );

    tmp<Field<Type> > ttoField
    (
        new Field<Type>(globalPatchA().globalPatch().size())
    );

    return ttoField;
}

template<class Type>
bool regionInterfaceType::couplesField(const word& fldName) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if
    (
        this->meshA().foundObject<fieldType>(fldName)
     && this->meshB().foundObject<fieldType>(fldName)
    )
    {
        return true;
    }

    return false;
}

template<class Type>
scalarField regionInterfaceType::rawJumpRes(const word& fldName) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    // Check if field of type fieldType and name fldName exists
    // on both sides. If not -> error
    if
    (
        !(this->couplesField<Type>(fldName))
    )
    {
        FatalErrorIn
            (
                "scalarField regionInterfaceType::rawJumpRes"
                "("
                    "const word& fldName"
                ") const"
            )
            << "Field " << fldName << " of type " << fieldType::typeName
            << " does not exist on both sides of the interface "
            << this->interfaceName()
            << abort(FatalError);
    }

    // If interface is not coupled return a zero residual
    if (!(this->coupled()))
    {
        scalarField zeroField(this->patchA().size(), 0.0);
        return zeroField;
    }

    const fvPatchField<Type>& patchAField =
            this->patchA().lookupPatchField<fieldType, Type>(fldName);

    const fvPatchField<Type>& patchBField =
            this->patchB().lookupPatchField<fieldType, Type>(fldName);

    if
    (
        isA<genericRegionCoupledJumpFvPatchField<Type>>(patchAField)
    )
    {
        const genericRegionCoupledJumpFvPatchField<Type>& jumpPatchAField =
            refCast<const genericRegionCoupledJumpFvPatchField<Type> >
            (patchAField);

        return jumpPatchAField.rawResidual();
    }
    else if
    (
        isA<genericRegionCoupledJumpFvPatchField<Type>>(patchBField)
    )
    {
        const genericRegionCoupledJumpFvPatchField<Type>& jumpPatchBField =
            refCast<const genericRegionCoupledJumpFvPatchField<Type> >
            (patchBField);

        const genericRegionCoupledFluxFvPatchField<Type>& fluxPatchAField =
            refCast<const genericRegionCoupledFluxFvPatchField<Type> >
            (patchAField);

        tmp<scalarField> resField =
            fluxPatchAField.interpolateFromNbrField
            (
                jumpPatchBField.rawResidual()
            );

        return resField();
    }

    FatalErrorIn
        (
            "scalarField regionInterfaceType::rawJumpRes"
            "("
                "const word& fldName"
            ") const"
        )
        << "Neither patch " << this->patchA().name()
        << "nor patch " << this->patchA().name()
        << "composing the interface " << this->interfaceName()
        << "are of type genericRegionCoupledJumpFvPatchField"
        << abort(FatalError);

    scalarField zeroField(this->patchA().size(), 0.0);
    return zeroField;
}


template<class Type>
scalarField regionInterfaceType::rawFluxRes(const word& fldName) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    // Check if field of type fieldType and name fldName exists
    // on both sides. If not -> error
    if
    (
        !(this->couplesField<Type>(fldName))
    )
    {
        FatalErrorIn
            (
                "scalarField regionInterfaceType::rawFluxRes"
                "("
                    "const word& fldName"
                ") const"
            )
            << "Field " << fldName << " of type " << fieldType::typeName
            << " does not exist on both sides of the interface "
            << this->interfaceName()
            << abort(FatalError);
    }

    // If interface is not coupled return a zero residual
    if (!(this->coupled()))
    {
        scalarField zeroField(this->patchA().size(), 0.0);
        return zeroField;
    }

    const fvPatchField<Type>& patchAField =
            this->patchA().lookupPatchField<fieldType, Type>(fldName);

    const fvPatchField<Type>& patchBField =
            this->patchB().lookupPatchField<fieldType, Type>(fldName);

    if
    (
        isA<genericRegionCoupledFluxFvPatchField<Type>>(patchAField)
    )
    {
        const genericRegionCoupledFluxFvPatchField<Type>& fluxPatchAField =
            refCast<const genericRegionCoupledFluxFvPatchField<Type> >
            (patchAField);

        return fluxPatchAField.rawResidual();
    }
    else if
    (
        isA<genericRegionCoupledFluxFvPatchField<Type>>(patchBField)
    )
    {
        const genericRegionCoupledFluxFvPatchField<Type>& fluxPatchBField =
            refCast<const genericRegionCoupledFluxFvPatchField<Type> >
            (patchBField);

        const genericRegionCoupledJumpFvPatchField<Type>& jumpPatchAField =
            refCast<const genericRegionCoupledJumpFvPatchField<Type> >
            (patchAField);


        tmp<scalarField> resField =
            jumpPatchAField.interpolateFromNbrField
            (
                fluxPatchBField.rawResidual()
            );

        return resField();
    }

    FatalErrorIn
    (
        "scalarField regionInterfaceType::rawFluxRes"
        "("
            "const word& fldName"
        ") const"
    )
    << "Neither patch " << this->patchA().name()
    << "nor patch " << this->patchA().name()
    << "composing the interface " << this->interfaceName()
    << "are of type genericRegionCoupledFluxFvPatchField"
    << abort(FatalError);

    scalarField zeroField(this->patchA().size(), 0.0);
    return zeroField;
}

template<class Type>
scalar regionInterfaceType::normJumpRes(const word& fldName) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    // Check if field of type fieldType and name fldName exists
    // on both sides. If not -> error
    if
    (
        !(this->couplesField<Type>(fldName))
    )
    {
        FatalErrorIn
            (
                "scalar regionInterfaceType::normJumpRes"
                "("
                    "const word& fldName"
                ") const"
            )
            << "Field " << fldName << " of type " << fieldType::typeName
            << " does not exist on both sides of the interface "
            << this->interfaceName()
            << abort(FatalError);
    }

    // If interface is not coupled return a zero residual
    if (!(this->coupled()))
    {
        return 0.0;
    }

    const fvPatchField<Type>& patchAField =
        this->patchA().lookupPatchField<fieldType, Type>(fldName);

    const fvPatchField<Type>& patchBField =
        this->patchB().lookupPatchField<fieldType, Type>(fldName);

    if
    (
        isA<genericRegionCoupledJumpFvPatchField<Type>>(patchAField)
    )
    {
        scalar patchJumpRes =
            refCast<const genericRegionCoupledJumpFvPatchField<Type> >
            (patchAField).normResidual();

        return patchJumpRes;

    }
    else if
    (
        isA<genericRegionCoupledJumpFvPatchField<Type>>(patchBField)
    )
    {
        scalar patchJumpRes =
            refCast<const genericRegionCoupledJumpFvPatchField<Type> >
            (patchBField).normResidual();

        return patchJumpRes;

    }

    Warning << "Coupled patchField of field" << fldName
    << " is not derived from genericRegionCoupledJumpFvPatchField" << nl
    << " No DNA interface residual control possible!" <<nl
    << "Setting globalMaxJumpRes to maximum value"
    << endl;

    return GREAT;
}


template<class Type>
scalar regionInterfaceType::normFluxRes(const word& fldName) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    // Check if field of type fieldType and name fldName exists
    // on both sides. If not -> error
    if
    (
        !(this->couplesField<Type>(fldName))
    )
    {
        FatalErrorIn
            (
                "scalar regionInterfaceType::normFluxRes"
                "("
                    "const word& fldName"
                ") const"
            )
            << "Field " << fldName << " of type " << fieldType::typeName
            << " does not exist on both sides of the interface "
            << this->interfaceName()
            << abort(FatalError);
    }


    // If interface is not coupled return a zero residual
    if (!(this->coupled()))
    {
        return 0.0;
    }

    const fvPatchField<Type>& patchAField =
        this->patchA().lookupPatchField<fieldType, Type>(fldName);

    const fvPatchField<Type>& patchBField =
        this->patchB().lookupPatchField<fieldType, Type>(fldName);

    if
    (
        isA<genericRegionCoupledFluxFvPatchField<Type>>(patchAField)
    )
    {
        scalar patchFluxRes =
            refCast<const genericRegionCoupledFluxFvPatchField<Type> >
            (patchAField).normResidual();

        return patchFluxRes;

    }
    else if
    (
        isA<genericRegionCoupledFluxFvPatchField<Type>>(patchBField)
    )
    {
        scalar patchFluxRes =
            refCast<const genericRegionCoupledFluxFvPatchField<Type> >
            (patchBField).normResidual();

        return patchFluxRes;

    }

    Warning << "Coupled patchFields of field" << fldName
    << " is not derived from genericRegionCoupledFluxFvPatchField" << nl
    << " No DNA interface residual control possible!" <<nl
    << "Setting globalMaxFluxRes to maximum value"
    << endl;

    return GREAT;
}

//-Specialisations
template<>
tmp<Field<scalar> > regionInterfaceType::interpolateFacesFromA
(
    const Field<scalar>& fromField
) const;

template<>
tmp<Field<vector> > regionInterfaceType::interpolateFacesFromA
(
    const Field<vector>& fromField
) const;

template<>
tmp<Field<scalar> > regionInterfaceType::interpolateFacesFromB
(
    const Field<scalar>& fromField
) const;

template<>
tmp<Field<vector> > regionInterfaceType::interpolateFacesFromB
(
    const Field<vector>& fromField
) const;

template<>
tmp<Field<scalar> > regionInterfaceType::interpolatePointsFromA
(
    const Field<scalar>& fromField
) const;

template<>
tmp<Field<vector> > regionInterfaceType::interpolatePointsFromA
(
    const Field<vector>& fromField
) const;

template<>
tmp<Field<scalar> > regionInterfaceType::interpolatePointsFromB
(
    const Field<scalar>& fromField
) const;

template<>
tmp<Field<vector> > regionInterfaceType::interpolatePointsFromB
(
    const Field<vector>& fromField
) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
