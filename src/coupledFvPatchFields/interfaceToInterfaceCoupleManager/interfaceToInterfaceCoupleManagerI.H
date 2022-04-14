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

inline const Foam::fvPatch& Foam::interfaceToInterfaceCoupleManager::refPatch() const
{
    return patch_;
};


inline const Foam::word& Foam::interfaceToInterfaceCoupleManager::neighbourRegionName() const
{
    return neighbourRegionName_;
};


inline const Foam::word& Foam::interfaceToInterfaceCoupleManager::neighbourPatchName() const
{
    return neighbourPatchName_;
};


inline const Foam::word& Foam::interfaceToInterfaceCoupleManager::neighbourFieldName() const
{
     return neighbourFieldName_;
};


inline const Foam::fvMesh& Foam::interfaceToInterfaceCoupleManager::neighbourRegion() const
{
    return localRegion_.objectRegistry::parent()
        .lookupObject<fvMesh>(neighbourRegionName_);
}


inline Foam::label Foam::interfaceToInterfaceCoupleManager::neighbourPatchID() const
{
    return neighbourRegion().boundaryMesh().findPatchID(neighbourPatchName_);
}


inline const Foam::fvPatch& Foam::interfaceToInterfaceCoupleManager::neighbourPatch() const
{
    return neighbourRegion().boundary()[neighbourPatchID()];
}


template<class Type>
inline const Foam::fvPatchField<Type>&
Foam::interfaceToInterfaceCoupleManager::neighbourPatchField() const
{
    return neighbourPatch().fvPatch::lookupPatchField
        <GeometricField<Type, fvPatchField, volMesh>, Type>
            (neighbourFieldName_);
}


inline Foam::word& Foam::interfaceToInterfaceCoupleManager::neighbourRegionName()
{
    return neighbourRegionName_;
};


inline Foam::word& Foam::interfaceToInterfaceCoupleManager::neighbourPatchName()
{
    return neighbourPatchName_;
};


template<class Type>
Foam::tmp<Field<Type> > 
Foam::interfaceToInterfaceCoupleManager::interpolateFromNbrField
(
    const Field<Type>& fromField
) const
{
    // B is neighbour
    if (refPatch().name() == rgInterface().patchA().name())
    {
        Field<Type> globalFromFld = 
            rgInterface().globalPatchB().patchFaceToGlobal(fromField);

        Field<Type> interpolatedGlobalFld =
            rgInterface().transferFacesFromB<Type>(globalFromFld);

        return
        (
            rgInterface().globalPatchB().globalFaceToPatch
            (
                interpolatedGlobalFld
            )
        );
    }

    // A is neighbour
    Field<Type> globalFromFld = 
        rgInterface().globalPatchA().patchFaceToGlobal(fromField);

    Field<Type> interpolatedGlobalFld = 
        rgInterface().transferFacesFromA<Type>(globalFromFld);

    return
    (
        rgInterface().globalPatchA().globalFaceToPatch
        (
            interpolatedGlobalFld
        )
    );
}


// ************************************************************************* //
