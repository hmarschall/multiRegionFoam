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

Class
    Foam::regionTypeTemplates

Description
    Template specialisations

SourceFiles
    regionTypeTemplates.C

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "fvMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class T>
fvMatrix<T>& regionType::getCoupledEqn
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
