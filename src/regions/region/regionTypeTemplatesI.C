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
    Template functions

SourceFiles
    regionTypeTemplatesI.C

\*---------------------------------------------------------------------------*/

#include "dynamicFvMesh.H"
#include "fvMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template <class T>
autoPtr<T> regionType::lookupOrRead
(
    const fvMesh& mesh,
    const word& fldName,
    const bool& read,
    const bool& write,
    const tmp<T> fld
)
{
    autoPtr<T> vfPtr(nullptr);

    if (mesh.foundObject<T>(fldName))
    {
        vfPtr.reset
        (
            const_cast<T*>
            (   
                &mesh.lookupObject<T>(fldName)
            )
        );
    }
    else
    {
        IOobject header
        (
            fldName,
            mesh.time().timeName(),
            mesh.time(),
            IOobject::NO_READ
        );

        IOobject io
        (
            fldName,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (!write)
        {
            io.writeOpt() = IOobject::NO_WRITE;
        }

        if
        (
            !read
         && !(header.headerOk())
        )
        {
            io.readOpt() = IOobject::NO_READ;

            vfPtr.reset
            (
                new T
                (
                    io,
                    fld
                )
            );
        }
        else // only if field exists and can be read (for restart)
        {
            vfPtr.reset
            (
                new T
                (
                    io,
                    mesh
                )
            );
        }
    }

    return vfPtr;
}

template <class T>
autoPtr<T> regionType::lookupOrRead
(
    const fvMesh& mesh,
    const word& fldName,
    const dimensioned<typename T::cmptType> dimVal,
    const bool& write
)
{
    autoPtr<T> vfPtr(nullptr);

    if (mesh.foundObject<T>(fldName))
    {
        vfPtr.reset
        (
            const_cast<T*>
            (   
                &mesh.lookupObject<T>(fldName)
            )
        );
    }
    else
    {
        IOobject io
        (
            fldName,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        );

        if (!write)
        {
            io.writeOpt() = IOobject::NO_WRITE;
        }

        vfPtr.reset
        (
            new T
            (
                io,
                mesh,
                dimVal
            )
        );
    }

    return vfPtr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
