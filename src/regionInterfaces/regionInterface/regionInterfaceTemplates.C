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

template<class Type>
tmp<Field<Type> > Foam::regionInterface::transferFacesFromA
(
    const Field<Type>& fromField
) const
{
//    Field<Type> toField(patch().size(), pTraits<Type>::zero);
    tmp<Field<Type> > ttoField(new Field<Type>(patch().size()));
    Field<Type>& toField = ttoField();

    transferFacesZoneToZone
    (
        meshA().name(),                     // from region name
        meshB().name(),                     // to region name
        globalPatchA().globalPatch(),       // from zone
        globalPatchB().globalPatch(),       // to zone
        fromField,                          // from field
        toField                             // to field
    );

//    return tmp<Field<Type> >(toField);
    return ttoField;
}

template<class Type>
void Foam::regionInterface::transferFacesAToB
(
    const Field<Type>& fromField,
    Field<Type>& toField
) const
{
    transferFacesZoneToZone
    (
        meshA().name(),                     // from region name
        meshB().name(),                     // to region name
        globalPatchA().globalPatch(),       // from zone
        globalPatchB().globalPatch(),       // to zone
        fromField,                          // from field
        toField                             // to field
    );
}

template<class Type>
void Foam::regionInterface::transferFacesBToA
(
    const Field<Type>& fromField,
    Field<Type>& toField
) const
{
    transferFacesZoneToZone
    (
        meshB().name(),                     // from region name
        meshA().name(),                     // to region name
        globalPatchB().globalPatch(),       // from zone
        globalPatchA().globalPatch(),       // to zone
        fromField,                          // from field
        toField                             // to field
    );
}

template<class Type>
void Foam::regionInterface::transferFacesZoneToZone
(
    const word& fromRegion,          // from region name
    const word& toRegion,            // to region name
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<Type>& fromField,    // from field
    Field<Type>& toField             // to field
) const
{
    // Check field sizes are correct
    if (fromField.size() != fromZone.size())
    {
        FatalErrorIn
        (
            "Foam::tmp< Field<Type> >\n"
            "Foam::regionInterface::transferFacesZoneToZone\n"
            "(\n"
            "    const standAlonePatch& fromZone,\n"
            "    const standAlonePatch& toZone,\n"
            "    const Field<Type>& fromField,\n"
            "    Field<Type>& toField\n"
            ") const"
        )   << "fromField is wrong size!" << nl
            << "fromField size: " << fromField.size()
            << ", fromZone size: " << fromZone.size()
            << abort(FatalError);
    }

    if (toField.size() != toZone.size())
    {
        FatalErrorIn
        (
            "Foam::tmp< Field<Type> >\n"
            "Foam::regionInterface::transferFacesZoneToZone\n"
            "(\n"
            "    const standAlonePatch& fromZone,\n"
            "    const standAlonePatch& toZone,\n"
            "    const Field<Type>& fromField,\n"
            "    Field<Type>& toField\n"
            ") const"
        )   << "toField is wrong size!" << nl
            << "toField size: " << toField.size()
            << ", toZone size: " << toZone.size()
            << abort(FatalError);
    }

    if (transferMethod_ == directMap)
    {
        FatalErrorIn
        (
            "Foam::tmp< Field<Type> >\n"
            "Foam::regionInterface::transferFacesZoneToZone\n"
            "(\n"
            "    const standAlonePatch& fromZone,\n"
            "    const standAlonePatch& toZone,\n"
            "    const Field<Type>& fromField,\n"
            "    Field<Type>& toField\n"
            ") const"
        )   << "Transfer method not implemented."
            << abort(FatalError);

//        if (fromRegion == fluidMesh().name() && toRegion == solidMesh().name())
//        {
//            const labelList& fluidToSolidMap = fluidToSolidFaceMap();
//            forAll(toField, faceI)
//            {
//                toField[faceI] = fromField[fluidToSolidMap[faceI]];
//            }
//        }
//        else if
//        (
//            fromRegion == solidMesh().name() && toRegion == fluidMesh().name()
//        )
//        {
//            const labelList& solidToFluidMap = solidToFluidFaceMap();
//            forAll(toField, faceI)
//            {
//                toField[faceI] = fromField[solidToFluidMap[faceI]];
//            }
//        }
//        else
//        {
//            FatalErrorIn
//            (
//                "Foam::tmp< Field<Type> >\n"
//                "Foam::regionInterface::transferFacesZoneToZone\n"
//                "(\n"
//                "    const standAlonePatch& fromZone,\n"
//                "    const standAlonePatch& toZone,\n"
//                "    const Field<Type>& fromField,\n"
//                "    Field<Type>& toField\n"
//                ") const"
//            )   << "Unknown regions:  " << fromRegion
//                << " and/or " << toRegion << abort(FatalError);
//        }
    }
    else if (transferMethod_ == RBF)
    {
        FatalErrorIn
        (
            "Foam::tmp< Field<Type> >\n"
            "Foam::regionInterface::transferFacesZoneToZone\n"
            "(\n"
            "    const standAlonePatch& fromZone,\n"
            "    const standAlonePatch& toZone,\n"
            "    const Field<Type>& fromField,\n"
            "    Field<Type>& toField\n"
            ") const"
        )   << "Transfer method not implemented."
            << abort(FatalError);

//        Info<< "Interpolating from " << fromRegion << " to " << toRegion
//            << " using RBF interpolation" << endl;

//        matrix fromRbfField(fromField.size(), int(pTraits<Type>::nComponents));
//        matrix toRbfField(toField.size(), int(pTraits<Type>::nComponents));

//        for(int cmptI = 0; cmptI < pTraits<Type>::nComponents; cmptI++)
//        {
//            const scalarField fromFieldCmptI = fromField.component(cmptI);

//            forAll(fromField, faceI)
//            {
//                fromRbfField(faceI, cmptI) = fromFieldCmptI[faceI];
//            }
//        }

//        if (fromRegion == fluidMesh().name() && toRegion == solidMesh().name())
//        {
//            rbfFluidToSolid()->interpolate(fromRbfField, toRbfField);
//        }
//        else if
//        (
//            fromRegion == solidMesh().name() && toRegion == fluidMesh().name()
//        )
//        {
//            rbfSolidToFluid()->interpolate(fromRbfField, toRbfField);
//        }
//        else
//        {
//            FatalErrorIn
//            (
//                "Foam::tmp< Field<Type> >\n"
//                "Foam::regionInterface::transferFacesZoneToZone\n"
//                "(\n"
//                "    const standAlonePatch& fromZone,\n"
//                "    const standAlonePatch& toZone,\n"
//                "    const Field<Type>& fromField,\n"
//                "    Field<Type>& toField\n"
//                ") const"
//            )   << "Unknown regions:  " << fromRegion
//                << " and/or " << toRegion << abort(FatalError);
//        }

//        for(int cmptI = 0; cmptI < pTraits<Type>::nComponents; cmptI++)
//        {
//            scalarField toFieldCmptI(toField.size(), 0.0);

//            forAll(toField, faceI)
//            {
//                toFieldCmptI[faceI] = toRbfField(faceI, cmptI);
//            }

//            toField.replace(cmptI, toFieldCmptI);
//        }
    }
    else if (transferMethod_ == GGI)
    {
        if (debug)
        {
            Info<< "Interpolating from " << fromRegion << " to " << toRegion
                << " using GGI interpolation" << endl;
        }

        if (fromRegion == meshB().name() && toRegion == meshA().name())
        {
            // fluid is the master; solid is the slave
            toField = ggiInterpolator().masterToSlave(fromField);
        }
        else if
        (
            fromRegion == meshA().name() && toRegion == meshB().name()
        )
        {
            toField = ggiInterpolator().slaveToMaster(fromField);
        }
        else
        {
            FatalErrorIn
            (
                "Foam::tmp< Field<Type> >\n"
                "Foam::regionInterface::transferFacesZoneToZone\n"
                "(\n"
                "    const standAlonePatch& fromZone,\n"
                "    const standAlonePatch& toZone,\n"
                "    const Field<Type>& fromField,\n"
                "    Field<Type>& toField\n"
                ") const"
            )   << "Unknown regions:  " << fromRegion
                << " and/or " << toRegion << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn
        (
            "Foam::tmp< Field<Type> >\n"
            "Foam::regionInterface::transferFacesZoneToZone\n"
            "(\n"
            "    const standAlonePatch& fromZone,\n"
            "    const standAlonePatch& toZone,\n"
            "    const Field<Type>& fromField,\n"
            "    Field<Type>& toField\n"
            ") const"
        )   << "Unknown transferMethod:  "
            << interfaceTransferMethodNames_[transferMethod_] << nl
            << "Available transfer methods are: "
            << interfaceTransferMethodNames_
            << abort(FatalError);
    }
}

