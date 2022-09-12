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
#include "dynamicFvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "fvMatrix.H"

#include "regionTypeList.H"
#include "regionInterfaceList.H"

#include "dnaControl.H"


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

        //- Reference to runTime
        const Time& runTime_;

        //- lists of regions
        autoPtr<regionTypeList> regions_;

        //- list of region interfaces
        autoPtr<regionInterfaceList> interfaces_;

        //- Lists of PtrLists of coupled finite volume fields
        //  per coupling method
        List<PtrList<volScalarField> > scalarFlds_;
        List<PtrList<volVectorField> > vectorFlds_;
        List<PtrList<volSymmTensorField> > symmTensorFlds_;
        List<PtrList<volTensorField> > tensorFlds_;
        List<PtrList<volVector4Field> > vector4Flds_;

        //- List of unique lists of coupled field names
        //  per coupling method
        List<hashedWordList> fldNames_;

    // Private member functions

        //- Solve pressure-velocity system using PIMPLE
        void solvePIMPLE();

        //- assemble and solve blockMatrix system
        template<class T>
        void assembleCoupledFields
        (
            List<PtrList<GeometricField<T, fvPatchField, volMesh> > >& flds,
            const List<hashedWordList>& fldNms
        ) const;

        template< template<class> class M, class T>
        void assembleAndSolveCoupledMatrix
        (
            PtrList<GeometricField<T, fvPatchField, volMesh> >& flds,
            word fldName
        ) const;

        template< template<class> class M, class T>
        void assembleAndSolveEqns
        (
            word fldName
        ) const;

protected:

        friend regionTypeList;

        friend regionType;

public:

    // Constructors

        //- Construct from dynamicFvMesh
        multiRegionSystem(const Time& runTime);


    //- Destructor
    virtual ~multiRegionSystem();


    //- Coupling options
    enum couplingMethods
    {
        partitioned,
        coupled
    };

    static const NamedEnum<couplingMethods, 2> couplingMethodsNames_;


    // Member Functions

        void preSolve();

        void solve();

        void setCoupledEqns();

        scalar getMinDeltaT();
        
        void postSolve();

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
