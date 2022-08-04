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

Description
    A template file for a coded refinementSelection

SourceFiles
    codedFieldBoundsRefinementTemplate.C

Author
    Mohammed Elwardi Fadeli (elwardifadeli@gmail.com)

\*---------------------------------------------------------------------------*/

#ifndef codedFieldBoundsRefinementTemplate_H
#define codedFieldBoundsRefinementTemplate_H

#include "refinementSelection.H"
#include "Switch.H"
#include "IOdictionary.H"
#include "volMesh.H"
#include "GeometricFields.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class ${typeName} Declaration
\*---------------------------------------------------------------------------*/

class ${typeName}
:
    public refinementSelection
{
    // Private data

        //- Name of the field to refine
        word fieldName_;

        //- Target field
        mutable volScalarField field_;

        //- Lower bound for the field
        mutable scalar lowerBound_;

        //- Upper bound for the field
        mutable scalar upperBound_;

        //- Whether to use cell-point-cell interpolation to smooth out the field
        //  before selection. Switched off by default
        Switch cellPointCellSmoothing_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ${typeName}(const ${typeName}&);

        //- Disallow default bitwise assignment
        void operator=(const ${typeName}&);


public:

    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;

    //- Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct from components
        ${typeName}(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~${typeName}();


    // Member Functions

        // Selection of refinement/unrefinement candidates

            //- Return transferable list of cells to refine
            virtual Xfer<labelList> refinementCellCandidates() const;

            //- Return transferable list of split points to unrefine
            virtual Xfer<labelList> unrefinementPointCandidates() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
