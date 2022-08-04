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

\*---------------------------------------------------------------------------*/

#include "refinementDetails.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(refinementDetails, 0);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::refinementDetails::readDict()
{
    // Read separate updates switch
    separateUpdates_ = false;
        //refinementDict_.lookupOrDefault<Switch>("separateUpdates", false);

    // Read maxRefinementLevel
    maxRefinementLevel_ = readLabel(refinementDict_.lookup("maxRefinementLevel"));

    // Read edgeBasedConsistency
    refinementDict_.lookup("edgeBasedConsistency") >> edgeBasedConsistency_ ;

    // Read nRefinementBufferLayers
    nRefinementBufferLayers_ = readLabel(refinementDict_.lookup("nRefinementBufferLayers"));

    // Read nUnrefinementBufferLayers
    nUnrefinementBufferLayers_ = readLabel(refinementDict_.lookup("nUnrefinementBufferLayers"));
}

void Foam::refinementDetails::setRefinementSettings
(
    enhancedRefinement& ref
) const
{

    // Gain access and modify some private members of refinement
    ref.setMaxRefinementLevel(maxRefinementLevel_);
    ref.setEdgeBasedConsistency(edgeBasedConsistency_);
    ref.setNRefinementBufferLayers(nRefinementBufferLayers_);
    ref.setNUnrefinementBufferLayers(nUnrefinementBufferLayers_);

    if (refinement::debug)
    Info<< "Changed refinement settings for "
        << refinementDict().name()
        << " to:" << nl
        << "{" << nl
        << tab << "maxRefinementLevel: " << ref.maxRefinementLevel() << nl
        << tab << "edgeBasedConsistency: " << ref.edgeBasedConsistency() << nl
        << tab << "nRefinementBufferLayers: " << ref.nRefinementBufferLayers() << nl
        << tab << "nUnrefinementBufferLayers: " << ref.nUnrefinementBufferLayers() << nl
        << "}" << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementDetails::refinementDetails()
:
    refinementDict_(nullptr),
    separateUpdates_(false),
    curTimeIndex_(-1),
    maxRefinementLevel_(1),
    edgeBasedConsistency_(false),
    nRefinementBufferLayers_(1),
    nUnrefinementBufferLayers_(1)
{
}


Foam::refinementDetails::refinementDetails
(
    const Time& time,
    const IOobject& io,
    const objectRegistry& db,
    const dictionary& refinementDict
)
:
    refinementDict_(refinementDict),
    separateUpdates_
    (
        false
        //refinementDict_.lookupOrDefault<Switch>("separateUpdates", false)
    ),
    curTimeIndex_(-1),
    maxRefinementLevel_(readLabel(refinementDict_.lookup("maxRefinementLevel"))),
    edgeBasedConsistency_
    (
        refinementDict_.lookupOrDefault<Switch>("edgeBasedConsistency", true)
    ),
    nRefinementBufferLayers_
    (
        readScalar(refinementDict_.lookup("nRefinementBufferLayers"))
    ),
    nUnrefinementBufferLayers_
    (
        readScalar(refinementDict_.lookup("nUnrefinementBufferLayers"))
    )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::refinementDetails::~refinementDetails()
{}

// ************************************************************************* //
