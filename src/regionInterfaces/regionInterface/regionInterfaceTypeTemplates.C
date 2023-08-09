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

#include "regionInterfaceType.H"

// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

namespace Foam
{


template<>
tmp<Field<scalar> > regionInterfaceType::interpolateFacesFromA
(
    const Field<scalar>& fromField
) const
{
    Field<scalar> globalFromFld =
        globalPatchA().patchFaceToGlobal(fromField);

    Field<scalar> globalToField(globalPatchB().globalPatch().size());

    interfaceToInterface().transferFacesZoneToZone
    (
        globalPatchA().globalPatch(),       // from zone
        globalPatchB().globalPatch(),       // to zone
        globalFromFld,                      // from field
        globalToField                       // to field
    );

    return globalPatchB().globalFaceToPatch
        (
            globalToField
        );
}

template<>
tmp<Field<vector> > regionInterfaceType::interpolateFacesFromA
(
    const Field<vector>& fromField
) const
{
    Field<vector> globalFromFld =
        globalPatchA().patchFaceToGlobal(fromField);

    Field<vector> globalToField(globalPatchB().globalPatch().size());

    interfaceToInterface().transferFacesZoneToZone
    (
        globalPatchA().globalPatch(),       // from zone
        globalPatchB().globalPatch(),       // to zone
        globalFromFld,                      // from field
        globalToField                       // to field
    );

    return globalPatchB().globalFaceToPatch
        (
            globalToField
        );
}



template<>
tmp<Field<scalar> > regionInterfaceType::interpolateFacesFromB
(
    const Field<scalar>& fromField
) const
{
    Field<scalar> globalFromFld =
        globalPatchB().patchFaceToGlobal(fromField);

    Field<scalar> globalToField(globalPatchA().globalPatch().size());

    interfaceToInterface().transferFacesZoneToZone
    (
        globalPatchB().globalPatch(),       // from zone
        globalPatchA().globalPatch(),       // to zone
        globalFromFld,                      // from field
        globalToField                       // to field
    );

    return globalPatchA().globalFaceToPatch
        (
            globalToField
        );
}

template<>
tmp<Field<vector> > regionInterfaceType::interpolateFacesFromB
(
    const Field<vector>& fromField
) const
{
    Field<vector> globalFromFld =
        globalPatchB().patchFaceToGlobal(fromField);

    Field<vector> globalToField(globalPatchA().globalPatch().size());

    interfaceToInterface().transferFacesZoneToZone
    (
        globalPatchB().globalPatch(),       // from zone
        globalPatchA().globalPatch(),       // to zone
        globalFromFld,                      // from field
        globalToField                       // to field
    );

    return globalPatchA().globalFaceToPatch
        (
            globalToField
        );
}

}

