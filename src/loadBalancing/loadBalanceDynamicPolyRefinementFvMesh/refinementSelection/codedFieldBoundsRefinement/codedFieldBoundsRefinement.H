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

Class
    Foam::codedFieldBoundsRefinement

Description
    Constructs on-the-fly refinement selections

    The hook functions take no arguments, should not return anything
    and get executed at the beginning of the real functions:

    // Executed when selecting refinement candidates
    void code()

    // (optional) Executed when selecting unrefinement candidates
    void localCode()

    but you have access to mutable class members :
        field_ is the target field created as fieldName
        lowerBound_
        uppderBound_
        mesh() gets const access to mesh object

Usage
    Example usage in dynamicMeshDict.refinements (you can use any number of
    these):
    \verbatim
    basedOnAlpha
    {
        refineInterval   1;
        unrefineInterval 1;
        maxRefinementLevel   3;
        separateUpdates false;

        refinementSelection
        {
            type        codedFieldBoundsRefinement;
            fieldName   target; // Name of the field to be created
            lowerBound  0; // Initial lower bound
            upperBound  1; // Initial upper bound

            cellPointCellSmoothing on;

            // Optional: include extra headers
            codeInclude
            #{
                #include <iostream>
            #}

            // Mandatory: Execute this before selecting refinement candidates
            code
            #{
                Info << "** dynamicCode **" << endl;
                lowerBound_ = 0.5;
                auto alpha1 = mesh().lookupObject<volScalarField>("alpha1");
                field_ == alpha1;
            #};

            // Optional: Execute this before selecting unrefinement candidates
            localCode
            #{
                Info << "** Rarely needed **" << endl;
            #}
        }
    }
    \endverbatim

SourceFiles
    codedFieldBoundsRefinement.C

Author
    Mohammed Elwardi Fadeli (elwardifadeli@gmail.com)

\*---------------------------------------------------------------------------*/

#ifndef codedFieldBoundsRefinement_H
#define codedFieldBoundsRefinement_H

#include "refinementSelection.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class codedFieldBoundsRefinement Declaration
\*---------------------------------------------------------------------------*/

class codedFieldBoundsRefinement
:
    public refinementSelection,
    public codedBase
{
protected:

    // Protected data

        //- Class name for generated type and library directory name
        word name_;

        //- Copy of the full dictionary containing the coeffs sub-dict
        dictionary dict_;

        //- Pointer to the real refinement selection object
        mutable autoPtr<refinementSelection> redirectRefinementSelectionPtr_;


    // Protected Member Functions

        //- Name of the class
        virtual const word& codeName() const
        {
            return name_;
        }

        //- Adapt the context for the current object
        virtual void prepare
        (
            dynamicCode&,
            const dynamicCodeContext&
        ) const;

        //- Return loaded libs
        virtual dlLibraryTable& libs() const;

        // Return a description (type + name) for the output
        virtual string description() const;

        // Get the dictionary to initialize the codeContext
        virtual const dictionary& codeDict() const;

        // Clear any redirected objects
        virtual void clearRedirect() const
        {
            redirectRefinementSelectionPtr_.clear();
        }

        //- Disallow default bitwise copy construct
        codedFieldBoundsRefinement(const codedFieldBoundsRefinement&);

        //- Disallow default bitwise assignment
        void operator=(const codedFieldBoundsRefinement&);

public:

    //- Static Data members

        //- Name of the C code template to be used
        static constexpr const char* const codeTemplateC
             = "codedFieldBoundsRefinementTemplate.C";

        //- Name of the H code template to be used
        static constexpr const char* const codeTemplateH
             = "codedFieldBoundsRefinementTemplate.H";

    //- Runtime type information
    TypeName("codedFieldBoundsRefinement");


    // Constructors

        //- Construct from components
        codedFieldBoundsRefinement(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~codedFieldBoundsRefinement();

    // Member Functions

        //- Compiled refinementSelection
        refinementSelection& redirectRefinementSelection() const;

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
