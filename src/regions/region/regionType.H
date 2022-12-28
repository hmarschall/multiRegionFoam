/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::regionType.C

SourceFiles
    regionType.C

\*---------------------------------------------------------------------------*/

#ifndef regionType_H
#define regionType_H

#include "IOdictionary.H"
#include "dimensionedScalar.H"

#include "Map.H"
#include "runTimeSelectionTables.H"

#include "volFields.H"
#include "dynamicFvMesh.H"
#include "fvMatrix.H"
#include "fvCFD.H"
#include "fvBlockMatrix.H"
#include "HashPtrTable.H"
#include "hashedWordList.H"

#include "mergePolyMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class multiRegionSystem;

/*---------------------------------------------------------------------------*\
                           Class regionType Declaration
\*---------------------------------------------------------------------------*/

class regionType
:
    public IOdictionary
{
    // Private data

        //- Mesh
        autoPtr<dynamicFvMesh> meshPtr_;

protected:

    // Protected data

        //- coupled matrix systems
        HashPtrTable<fvMatrix<scalar> > fvScalarMatrices;
        HashPtrTable<fvMatrix<vector> > fvVectorMatrices;
        HashPtrTable<fvMatrix<symmTensor> > fvSymmTensorMatrices;
        HashPtrTable<fvMatrix<tensor> > fvTensorMatrices;
        HashPtrTable<fvBlockMatrix<vector4> > fvVector4Matrices;

    // Protected member functions

        //- Lookup or read field
        template <class T>
        inline autoPtr<T> lookupOrRead
        (
            const fvMesh& mesh,
            const word& fldName,
            const bool& read=true,
            const bool& write=true,
            const tmp<T> fld = tmp<T>(nullptr)
        );

        template <class T>
        inline autoPtr<T> lookupOrRead
        (
            const fvMesh& mesh,
            const word& fldName,
            const dimensioned<typename T::cmptType> dimVal,
            const bool& write=true
        );

public:

    //- Runtime type information
    TypeName("regionType");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            regionType,
            dictionary,
            (
                const Time& runTime,
                const word& name
            ),
            (runTime, name)
        );

    // Constructors
        regionType
        (
            const Time& runTime,
            const word& name
        );

    //- Destructor
    virtual ~regionType(){}

        // Selectors
        static autoPtr<regionType> New
        (
            const Time& runTime,
            const word& name,
            const word& modelType
        );

    // Member Functions

        //- Return reference to the mesh
        dynamicFvMesh& mesh()
        {
            return meshPtr_();
        }

        //- Return const reference to the mesh
        const fvMesh& mesh() const
        {
            return meshPtr_();
        }

        //- Return time
        const Time& runTime() const
        {
            return meshPtr_().time();
        }

        bool update()
        {
            return mesh().update();
        }

        // virtual functions

        virtual word regionTypeName() const
        {
            return word("regionTypeBase");
        }

        virtual bool active() const
        {
            return true;
        }

        virtual void correct() = 0;

        virtual scalar getMinDeltaT() = 0;

        virtual void solveRegion() = 0;

        virtual void prePredictor() = 0;

        virtual void momentumPredictor() = 0;

        virtual void pressureCorrector() = 0;

        virtual void meshMotionCorrector() = 0;

        virtual void setCoupledEqns() = 0;
        
        virtual void postSolve() = 0;

        template<class T>
        void relaxEqn(fvMatrix<T>& eqn)
        {
            eqn.relax();
        }

        template<class T>
        void relaxEqn(fvBlockMatrix<T>& eqn)
        {
            //TODO: Relaxation of fvBlockMatricies currently not supported
        }

        //virtual void clear() = 0;

        // Get a ref to the fvMatrix 
        template< template<class> class M, class T>
        M<T>& getCoupledEqn
        (
            word name
        );

        // Get a ref to the fvMatrix 
        template<class T>
        bool clearCoupledEqn
        (
            const T& fld,
            const word& typeName
        );

        // Get a ref to the fields 
        template<class T>
        const T& getObject(word name, const fvMesh& mesh);

        // Check if region holds fvMatrix 
        bool foundCoupledEqn
        (
            word name
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "regionTypeTemplatesI.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
