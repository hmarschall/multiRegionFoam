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
#include "fvMesh.H"
#include "fvMatrix.H"
#include "HashPtrTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class multiRegionSystem;

/*---------------------------------------------------------------------------*\
                           Class regionType Declaration
\*---------------------------------------------------------------------------*/

class regionType
:
    public fvMesh
{
public:

    typedef
        HashTable
        <
            Foam::tmp<Foam::fvScalarMatrix>, word, word::hash
        >
        coupledScalarEqns;

protected:

    // Protected data

        //- mesh
        const fvMesh& mesh_;

        const IOdictionary dict_;

        //- coupled governing equations
        HashPtrTable<fvMatrix<scalar> > fvScalarMatrices;
        HashPtrTable<fvMatrix<vector> > fvVectorMatrices;
        HashPtrTable<fvMatrix<symmTensor> > fvSymmTensorMatrices;
        HashPtrTable<fvMatrix<tensor> > fvTensorMatrices;

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
                const fvMesh& mesh,
                const word& name
            ),
            (mesh, name)
        );

    // Constructors
        regionType
        (
            const fvMesh& mesh,
            const word& name
        );

    //- Destructor
    virtual ~regionType();

        // Selectors
        static autoPtr<regionType> New
        (
            const fvMesh& mesh,
            const word& name,
            const word& modelType
        );

    // Member Functions

        // virtual functions

        virtual bool active() const
        {
            return true;
        }

        // virtual coupledScalarEqns coupledScalarEqnsTable() const;

        virtual void correct() = 0;

        virtual void setRDeltaT() = 0;

        virtual void solveRegion() = 0;

        virtual void setCoupledEqns() = 0;

        // Get a ref to the fvMatrix 
        template<class T>
        fvMatrix<T>& getCoupledEqn
        (
            word name
        );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
