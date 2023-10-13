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

#include "interfaceToInterfaceCoupleManager.H"
#include "OFstream.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::interfaceToInterfaceCoupleManager::assembleName
(
    const fvPatch& refPatch,
    const word& nbrRgName,
    const word& nbrPatchName,
    const word& typeName,
    const bool& reverseOrder
)
{
    const fvMesh& mesh = refPatch.boundaryMesh().mesh();

    word meshAName = nbrRgName;
    word patchAName = nbrPatchName;
    word meshBName = mesh.name();
    word patchBName = refPatch.name();

    if (reverseOrder)
    {
        meshAName = mesh.name();
        patchAName = refPatch.name();
        meshBName = nbrRgName;
        patchBName = nbrPatchName;
    }

    word PatchAName = word(toupper(patchAName[0]));
    PatchAName = PatchAName + patchAName.substr(1);

    word MeshBName = word(toupper(meshBName[0]));
    MeshBName = MeshBName + meshBName.substr(1);

    word PatchBName = word(toupper(patchBName[0]));
    PatchBName = PatchBName + patchBName.substr(1);

    word InterfaceTypeName = word(toupper(typeName[0]));
    InterfaceTypeName = InterfaceTypeName + typeName.substr(1);

    return
    (
        meshAName + PatchAName
      + MeshBName + PatchBName
      + InterfaceTypeName
    );
}

const Foam::regionInterfaceType&
Foam::interfaceToInterfaceCoupleManager::rgInterface() const
{
    const fvMesh& mesh = refPatch().boundaryMesh().mesh();
    const objectRegistry& obr = mesh.objectRegistry::parent();

    word rgIntName = assembleName
    (
        refPatch(),
        neighbourRegionName_,
        neighbourPatchName_,
        typeName_
    );

    word rgIntNameRev = assembleName
    (
        refPatch(),
        neighbourRegionName_,
        neighbourPatchName_,
        typeName_,
        true //reverse order
    );

    if
    (
        !obr.foundObject<regionInterfaceType>(rgIntName)
     && !obr.foundObject<regionInterfaceType>(rgIntNameRev)
    )
    {
        FatalErrorIn("interfaceToInterfaceCoupleManager::rgInterface()")
            << "regionInterfaceType object not found but required."
            << rgIntName << " " << rgIntNameRev
            << abort(FatalError);
    } else if
    (
        obr.foundObject<regionInterfaceType>(rgIntName)
     && obr.foundObject<regionInterfaceType>(rgIntNameRev)
    )
    {
        FatalErrorIn("interfaceCoupledVelocityValue::")
            << "regionInterfaceType object names ambiguous:"
            << rgIntName << " vs. " << rgIntNameRev << nl
            << "Choose unique patch/region names."
            << abort(FatalError);
    }

    if (obr.foundObject<regionInterfaceType>(rgIntNameRev))
    {
        return obr.lookupObject<regionInterfaceType>(rgIntNameRev);
    }

    return obr.lookupObject<regionInterfaceType>(rgIntName);
}

const Foam::fvMesh&
Foam::interfaceToInterfaceCoupleManager::nbrMesh() const
{
    if (refPatch().name() == rgInterface().patchA().name())
    {
        return rgInterface().meshB();
    }

    return rgInterface().meshA();
}

const Foam::fvPatch&
Foam::interfaceToInterfaceCoupleManager::nbrPatch() const
{
    if (refPatch().name() == rgInterface().patchA().name())
    {
        return rgInterface().patchB();
    }

    return rgInterface().patchA();
}

const Foam::Switch&
Foam::interfaceToInterfaceCoupleManager::interfaceCoupled() const
{
    return rgInterface().coupled();
}

void Foam::interfaceToInterfaceCoupleManager::updateRegionInterface()
{
    regionInterfaceType& rgInt = const_cast<regionInterfaceType&>
        (
            rgInterface()
        );

    rgInt.update();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceToInterfaceCoupleManager::interfaceToInterfaceCoupleManager
(
    const fvPatch& patch,
    const word type
)
:
    patch_(patch),
    typeName_(type),
    neighbourRegionName_(),
    neighbourPatchName_(),
    neighbourFieldName_(),
    localRegion_(patch_.boundaryMesh().mesh())
{}


Foam::interfaceToInterfaceCoupleManager::interfaceToInterfaceCoupleManager
(
    const fvPatch& patch,
    const dictionary& dict,
    const word type
)
:
    patch_(patch),
    typeName_(dict.lookup("interfaceType")),
    neighbourRegionName_(dict.lookup("neighbourRegionName")),
    neighbourPatchName_(dict.lookup("neighbourPatchName")),
    neighbourFieldName_(dict.lookup("neighbourFieldName")),
    localRegion_(patch_.boundaryMesh().mesh())
{}


Foam::interfaceToInterfaceCoupleManager::interfaceToInterfaceCoupleManager
(
    const interfaceToInterfaceCoupleManager& pcm
)
:
    patch_(pcm.refPatch()),
    typeName_(pcm.typeName_),
    neighbourRegionName_(pcm.neighbourRegionName()),
    neighbourPatchName_(pcm.neighbourPatchName()),
    neighbourFieldName_(pcm.neighbourFieldName()),
    localRegion_(patch_.boundaryMesh().mesh())
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceToInterfaceCoupleManager::~interfaceToInterfaceCoupleManager()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::interfaceToInterfaceCoupleManager::writeEntries(Ostream& os) const
{
    os.writeKeyword("neighbourRegionName");
    os << neighbourRegionName_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourPatchName");
    os << neighbourPatchName_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourFieldName");
    os << neighbourFieldName_ << token::END_STATEMENT << nl;
}

// ************************************************************************* //
