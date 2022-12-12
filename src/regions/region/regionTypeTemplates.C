/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

Description
    Template specialisations

SourceFiles
    regionTypeTemplates.C

\*---------------------------------------------------------------------------*/

#include "dynamicFvMesh.H"
#include "fvMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template< template<class> class M, class T>
M<T>& regionType::getCoupledEqn
(
    word name
)
{
    notImplemented
    (
        "regionTypeTemplates.C\n"
        "fvMatrix<T>& regionType::getCoupledEqn\n"
        "(\n"
        "word name\n"
        ")\n"
        "not implemented"
    );
}

template<>
fvMatrix<scalar>& regionType::getCoupledEqn
(
    word name
)
{
//    HashPtrTable<fvScalarMatrix>::iterator it = fvScalarMatrices.find(name);

//    if (it == fvScalarMatrices.end()) // not found
//    {
//        FatalErrorIn("regionType::getCoupledEqn")
//           << "Equation of type fvScalarMatrix "
//           << "with name " << name
//           << " not found"
//           << exit(FatalError);
//    }

//    return **it;

    return *fvScalarMatrices[name];
}

template<>
fvMatrix<vector>& regionType::getCoupledEqn
(
    word name
)
{
    return *fvVectorMatrices[name];
}

template<>
fvMatrix<symmTensor>& regionType::getCoupledEqn
(
    word name
)
{
    return *fvSymmTensorMatrices[name];
}

template<>
fvMatrix<tensor>& regionType::getCoupledEqn
(
    word name
)
{
    return *fvTensorMatrices[name];
}

template<>
fvBlockMatrix<vector4>& regionType::getCoupledEqn
(
    word name
)
{
    return *fvVector4Matrices[name];
}

bool regionType::foundCoupledEqn
(
    word name
)
{
    return 
    (
        fvScalarMatrices.found(name) ||
        fvVectorMatrices.found(name) ||
        fvSymmTensorMatrices.found(name) ||
        fvTensorMatrices.found(name) ||
        fvVector4Matrices.found(name)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
