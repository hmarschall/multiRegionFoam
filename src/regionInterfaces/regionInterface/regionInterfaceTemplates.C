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

#include "regionInterface.H"

// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

namespace Foam
{

template<class Type>
tmp<Field<Type> > regionInterface::transferFacesFromA
(
    const Field<Type>& fromField
) const
{
    notImplemented
    (
        "regionInterfaceTemplates.C\n"
        "tmp<Field<Type> > regionInterface::transferFacesFromA\n"
        "(\n"
        "const Field<Type>& fromField\n"
        ") const\n"
        "not implemented"
    );

    tmp<Field<Type> > ttoField(new Field<Type>(fromField.size()));

    return ttoField;
}

template<>
tmp<Field<scalar> > regionInterface::transferFacesFromA
(
    const Field<scalar>& fromField
) const
{
    tmp<Field<scalar> > ttoField(new Field<scalar>(fromField));
    Field<scalar>& toField = ttoField();

    interfaceToInterface().transferFacesZoneToZone
    (
        globalPatchA().globalPatch(),       // from zone
        globalPatchB().globalPatch(),       // to zone
        fromField,                          // from field
        toField                             // to field
    );

    return ttoField;
}

template<>
tmp<Field<vector> > regionInterface::transferFacesFromA
(
    const Field<vector>& fromField
) const
{
    tmp<Field<vector> > ttoField(new Field<vector>(fromField));
    Field<vector>& toField = ttoField();

    interfaceToInterface().transferFacesZoneToZone
    (
        globalPatchA().globalPatch(),       // from zone
        globalPatchB().globalPatch(),       // to zone
        fromField,                          // from field
        toField                             // to field
    );

    return ttoField;
}

template<class Type>
tmp<Field<Type> > regionInterface::transferFacesFromB
(
    const Field<Type>& fromField
) const
{
    notImplemented
    (
        "regionInterfaceTemplates.C\n"
        "tmp<Field<Type> > regionInterface::transferFacesFromB\n"
        "(\n"
        "const Field<Type>& fromField\n"
        ") const\n"
        "not implemented"
    );

    tmp<Field<Type> > ttoField(new Field<Type>(fromField.size()));

    return ttoField;
}

template<>
tmp<Field<scalar> > regionInterface::transferFacesFromB
(
    const Field<scalar>& fromField
) const
{
    tmp<Field<scalar> > ttoField(new Field<scalar>(fromField.size()));
    Field<scalar>& toField = ttoField();

    interfaceToInterface().transferFacesZoneToZone
    (
        globalPatchB().globalPatch(),       // from zone
        globalPatchA().globalPatch(),       // to zone
        fromField,                          // from field
        toField                             // to field
    );

    return ttoField;
}

template<>
tmp<Field<vector> > regionInterface::transferFacesFromB
(
    const Field<vector>& fromField
) const
{
    tmp<Field<vector> > ttoField(new Field<vector>(fromField.size()));
    Field<vector>& toField = ttoField();

    interfaceToInterface().transferFacesZoneToZone
    (
        globalPatchB().globalPatch(),       // from zone
        globalPatchA().globalPatch(),       // to zone
        fromField,                          // from field
        toField                             // to field
    );

    return ttoField;
}

}

