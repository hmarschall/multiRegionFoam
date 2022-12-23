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

// Get coupled equations from tables
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
    HashPtrTable<fvScalarMatrix>::iterator it =
        fvScalarMatrices.find(name);

    if (it == fvScalarMatrices.end()) // not found
    {
        FatalErrorIn("regionType::getCoupledEqn")
           << "Equation of type fvScalarMatrix "
           << "with name " << name
           << " not found"
           << exit(FatalError);
    }

    return **it;

//    return *fvScalarMatrices[name];
}

template<>
fvMatrix<vector>& regionType::getCoupledEqn
(
    word name
)
{
    HashPtrTable<fvVectorMatrix>::iterator it =
        fvVectorMatrices.find(name);

    if (it == fvVectorMatrices.end()) // not found
    {
        FatalErrorIn("regionType::getCoupledEqn")
           << "Equation of type fvVectorMatrix "
           << "with name " << name
           << " not found"
           << exit(FatalError);
    }

    return **it;

//    return *fvVectorMatrices[name];
}

template<>
fvMatrix<symmTensor>& regionType::getCoupledEqn
(
    word name
)
{
    HashPtrTable<fvSymmTensorMatrix>::iterator it =
        fvSymmTensorMatrices.find(name);

    if (it == fvSymmTensorMatrices.end()) // not found
    {
        FatalErrorIn("regionType::getCoupledEqn")
           << "Equation of type fvSymmTensorMatrix "
           << "with name " << name
           << " not found"
           << exit(FatalError);
    }

    return **it;

//    return *fvSymmTensorMatrices[name];
}

template<>
fvMatrix<tensor>& regionType::getCoupledEqn
(
    word name
)
{
    HashPtrTable<fvTensorMatrix>::iterator it =
        fvTensorMatrices.find(name);

    if (it == fvTensorMatrices.end()) // not found
    {
        FatalErrorIn("regionType::getCoupledEqn")
           << "Equation of type fvTensorMatrix "
           << "with name " << name
           << " not found"
           << exit(FatalError);
    }

    return **it;

//    return *fvTensorMatrices[name];
}

template<>
fvBlockMatrix<vector4>& regionType::getCoupledEqn
(
    word name
)
{
    HashPtrTable<fvBlockMatrix<vector4> >::iterator it =
        fvVector4Matrices.find(name);

    if (it == fvVector4Matrices.end()) // not found
    {
        FatalErrorIn("regionType::getCoupledEqn")
           << "Equation of type fvBlockMatrix<vector4> "
           << "with name " << name
           << " not found"
           << exit(FatalError);
    }

    return **it;

//    return *fvVector4Matrices[name];
}

// Clear coupled equations in tables
template<class T>
bool regionType::clearCoupledEqn
(
    const T& fld
)
{
    notImplemented
    (
        "regionTypeTemplates.C\n"
        "fvMatrix<T>& regionType::clearCoupledEqn\n"
        "(\n"
        "word name\n"
        ")\n"
        "not implemented"
    );

    return false;
}

template<>
bool regionType::clearCoupledEqn
(
    const volScalarField& fld
)
{
    HashPtrTable<fvScalarMatrix>::iterator it =
        fvScalarMatrices.find
        (
            fld.name() + fld.mesh().name() + "Eqn"
        );

    return fvScalarMatrices.erase(it);
}

template<>
bool regionType::clearCoupledEqn
(
    const volVectorField& fld
)
{
    HashPtrTable<fvVectorMatrix>::iterator it =
        fvVectorMatrices.find
        (
            fld.name() + fld.mesh().name() + "Eqn"
        );

    return fvVectorMatrices.erase(it);
}

template<>
bool regionType::clearCoupledEqn
(
    const volSymmTensorField& fld
)
{
    HashPtrTable<fvSymmTensorMatrix>::iterator it =
        fvSymmTensorMatrices.find
        (
            fld.name() + fld.mesh().name() + "Eqn"
        );

    return fvSymmTensorMatrices.erase(it);
}

template<>
bool regionType::clearCoupledEqn
(
    const volTensorField& fld
)
{
    HashPtrTable<fvTensorMatrix>::iterator it =
        fvTensorMatrices.find
        (
            fld.name() + fld.mesh().name() + "Eqn"
        );

    return fvTensorMatrices.erase(it);
}

template<>
bool regionType::clearCoupledEqn
(
    const volVector4Field& fld
)
{
    HashPtrTable<fvBlockMatrix<vector4> >::iterator it =
        fvVector4Matrices.find
        (
            fld.name() + fld.mesh().name() + "Eqn"
        );

    return fvVector4Matrices.erase(it);
}

// Find coupled equations in tables
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

// Get ref to registered objects by name
template<class T>
const T& getObject
(
    word name,
    const fvMesh& mesh
)
{
    notImplemented
    (
        "regionTypeTemplates.C\n"
        "fvMatrix<T>& regionType::getObject\n"
        "(\n"
        "const word& name, const fvMesh& mesh\n"
        ")\n"
        "not implemented"
    );

    return
    (
        mesh.thisDb().lookupObject<T>(name)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
