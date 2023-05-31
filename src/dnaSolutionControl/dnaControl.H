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
    Foam::dnaControl

Description
    Class for controling the Dirichlet-Neumann iteration procedure based on
    the interface residuals

\*---------------------------------------------------------------------------*/

#ifndef dnaControl_H
#define dnaControl_H


#include "IOdictionary.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "HashPtrTable.H"
#include "regionInterfaceTypeList.H"
#include "genericRegionCoupledFluxFvPatchField.H"
#include "genericRegionCoupledJumpFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class dnaControl Declaration
\*---------------------------------------------------------------------------*/

class dnaControl
:
    public IOobject
{

public:

    struct fieldData
    {
        word name;
        scalar maxJumpRes;
        scalar maxFluxRes;
        Switch outputJumpResField;
        Switch outputFluxResField;
    };

protected:

    // Protected data

        //-Field name for which dna is executed
        word fldName_;

        //- multiRegion properties dictionary
        //- DNA solution 
        IOdictionary multiRegionProperties_;

        //- reference to list of region interfaces
        const regionInterfaceTypeList& interfaces_;

        //- List of residual data per sub field
        List<fieldData> dnaResControl_;

        //- List of max interface residuals for each subField
        List<scalar> globalMaxJumpRes_;
        List<scalar> globalMaxFluxRes_;

        //- Hashed ptr tables for residual fields
        HashPtrTable<volScalarField> initJumpResFlds_;
        HashPtrTable<volScalarField> initFluxResFlds_;
        HashPtrTable<volScalarField> finalJumpResFlds_;
        HashPtrTable<volScalarField> finalFluxResFlds_;

        //- Maximum coupling iterations (partitioned)
        label maxCoupleIter_;

        //- Current corrector loop index
        label corr_;


    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        virtual void read();

        template<class Type>
        void maxTypeRes
        (
            const regionInterfaceType& interface,
            const word& fldName,
            scalar& globalMaxJumpRes,
            scalar& globalMaxFluxRes
        );

        void maxRes
        (
            const regionInterfaceType& interface,
            const word& fldName,
            scalar& globalMaxJumpRes,
            scalar& globalMaxFluxRes
        );

        template<class Type>
        void writeResFlds
        (
            const regionInterfaceType& interface,
            const word& fldName,
            const Switch outputJumpResField,
            const Switch outputFluxResField,
            bool final
        );

        void writeResFlds
        (
            const regionInterfaceType& interface,
            const word& fldName,
            const Switch outputJumpResField,
            const Switch outputFluxResField,
            bool final
        );

        void outputMaxResInfo
        (
            bool criteriaSatisfied,
            label corr
        );

        //- Return true if all convergence checks are satisfied
        virtual bool criteriaSatisfied();

private:

    // Private member functions

        //- Disallow default bitwise copy construct
        dnaControl(const dnaControl&);

        //- Disallow default bitwise assignment
        void operator=(const dnaControl&);

        //- create residual fields that should be written to disk
        void createResFlds();

public:

    // Static Data Members

        //- Run-time type information
        TypeName("dnaControl");


    // Constructors

        //- Construct from mesh and the name of control sub-dictionary
        dnaControl
        (
            const Time& runTime,
            const word& fldName,
            const regionInterfaceTypeList& interfaces
        );


    //- Destructor
    virtual ~dnaControl();

    // Member Functions

        // Evolution

            //- DNA loop
            virtual bool loop();
    
        // Access

            //- Return the solution dictionary
            inline const dictionary& dict() const;

            //- Current corrector loop index
            inline label corr() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dnaControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
//     #include "solutionControlTemplates.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //