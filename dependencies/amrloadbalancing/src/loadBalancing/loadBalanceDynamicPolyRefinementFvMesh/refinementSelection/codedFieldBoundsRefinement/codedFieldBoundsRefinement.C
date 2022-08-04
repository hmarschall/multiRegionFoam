/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    Mohammed Elwardi Fadeli (elwardifadeli@gmail.com)

\*---------------------------------------------------------------------------*/

#include "codedFieldBoundsRefinement.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "volPointInterpolation.H"
#include "pointVolInterpolation.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(codedFieldBoundsRefinement, 0);
addToRunTimeSelectionTable
(
    refinementSelection,
    codedFieldBoundsRefinement,
    dictionary
);

}

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

void Foam::codedFieldBoundsRefinement::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    // Filtering vars
    dynCode.setFilterVariable("typeName", name_);

    // Compile filtered file and add header
    dynCode.addCompileFile(codeTemplateC);
    dynCode.addCopyFile(codeTemplateH);

    // Make/options 
    dynCode.setMakeOptions
    (
        "EXE_INC = -g \\\n"
        "-I$(LIB_SRC)/finiteVolume/lnInclude \\\n"
        "-I$(LIB_SRC)/dynamicMesh/dynamicMesh/lnInclude \\\n"
        "-I$(LIB_SRC)/dynamicMesh/dynamicFvMesh/lnInclude \\\n"
        "-I$(LIB_SRC)/dynamicMesh/topoChangerFvMesh/lnInclude \\\n"
        "-I$(LIB_SRC)/dynamicMesh/dynamicTopoFvMesh/lnInclude \\\n"
      + context.options()
      + "\n\nLIB_LIBS = \\\n"
        "    -lfiniteVolume \\\n"
        "    -ldynamicMesh \\\n"
        "    -ldynamicFvMesh \\\n"
        "    -ltopoChangerFvMesh \\\n"
      + context.libs()
    );
}

Foam::dlLibraryTable& Foam::codedFieldBoundsRefinement::libs() const
{
    return const_cast<Foam::dlLibraryTable&>(mesh().time().libs());
}

Foam::string Foam::codedFieldBoundsRefinement::description() const
{
    return string("refinementSelection::") + name_;
}

const Foam::dictionary& Foam::codedFieldBoundsRefinement::codeDict() const
{
    return coeffDict();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::codedFieldBoundsRefinement::codedFieldBoundsRefinement
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    refinementSelection(mesh, dict),
    codedBase(),
    name_(),
    dict_(dict),
    redirectRefinementSelectionPtr_()
{
    // Do this here once, other calls to redirectRefinementSelection
    // will only check for the pointer's validity
    word dictName = dict_.dictName();
    name_ = "codedRefinement_" + dictName.replaceAll(':', "");
    updateLibrary(name_);
    redirectRefinementSelection();
}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::codedFieldBoundsRefinement::~codedFieldBoundsRefinement()
{}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

Foam::refinementSelection&
Foam::codedFieldBoundsRefinement::redirectRefinementSelection() const
{
    if (!redirectRefinementSelectionPtr_.valid())
    {
        dictionary constructDict(dict_);
        constructDict.subDict("refinementSelection").set("type", name_);
        redirectRefinementSelectionPtr_.reset
        (
            refinementSelection::New
            (
                mesh(),
                constructDict
            ).ptr()
        );
    }
    return redirectRefinementSelectionPtr_();
}


Foam::Xfer<Foam::labelList>
Foam::codedFieldBoundsRefinement::refinementCellCandidates() const
{
    if (debug)
    {
        Info<< "codedFieldBoundsRefinement::refinementCellCandidates"
            << " for selection " << name_ << endl;
    }

    updateLibrary(name_);
    return redirectRefinementSelection().refinementCellCandidates();
}


Foam::Xfer<Foam::labelList>
Foam::codedFieldBoundsRefinement::unrefinementPointCandidates() const
{
    if (debug)
    {
        Info<< "codedFieldBoundsRefinement::unrefinementPointCandidates"
            << " for selection " << name_ << endl;
    }

    updateLibrary(name_);
    return redirectRefinementSelection().unrefinementPointCandidates();
}


// ************************************************************************* //
