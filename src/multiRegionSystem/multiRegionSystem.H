/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::multiRegionSystem

Description

SourceFiles
    multiRegionSystem.C

\*---------------------------------------------------------------------------*/

#ifndef multiRegionSystem_H
#define multiRegionSystem_H

#include "IOReferencer.H"
#include "IOdictionary.H"
#include "hashedWordList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "fvMatrix.H"

#include "regionTypeList.H"
#include "regionInterfaceList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class regionTypeList;
class regionInterfaceList;

/*---------------------------------------------------------------------------*\
                      Class multiRegionSystem Declaration
\*---------------------------------------------------------------------------*/

class multiRegionSystem
:
    public IOdictionary
{
    // Private data

        typedef
            HashTable
            <
                wordList, interfaceKey, interfaceKey::hash
            >
            fieldsTable;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- lists of regions
        autoPtr<regionTypeList> regions_;

        //- list of region interfaces
        autoPtr<regionInterfaceList> interfaces_;

        //- PtrLists of finite volume fields
        PtrList<volScalarField> scalarFlds_;
        PtrList<volVectorField> vectorFlds_;
        PtrList<volSymmTensorField> symmTensorFlds_;
        PtrList<volTensorField> tensorFlds_;

        //- Maximum coupling iterations (partitioned)
        label maxCoupleIter_;

    // Private member functions

        //- assemble and solve blockMatrix system
        template<class T>
        void assembleCoupledFields
        (
            PtrList<GeometricField<T, fvPatchField, volMesh> >& flds,
            regionTypeList& regions,
            const hashedWordList& fldNames
        ) const;

        template<class T>
        void assembleAndSolveBlockMatrix
        (
            PtrList<GeometricField<T, fvPatchField, volMesh> >& flds,
            word fldName
        ) const;

        template<class T>
        void assembleAndSolveEqns
        (
            PtrList<GeometricField<T, fvPatchField, volMesh> >& flds,
            regionTypeList& regions,
            word fldName
        ) const;

protected:

        friend regionTypeList;

        friend regionType;


public:

    // Constructors

        //- Construct from fvMesh
        multiRegionSystem(const fvMesh&);


    //- Destructor
    virtual ~multiRegionSystem();


    // Member Functions

        // Correct
        void correct();

        void solve();

        void setCoupledEqns();

        void setRDeltaT();

        // Access

        regionTypeList& regions();

        regionInterfaceList& interfaces();

        //- Read base dictionary
//        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
