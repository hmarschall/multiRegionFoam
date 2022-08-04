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

Description
    A set of macro to make writing tests easier

Author
    Mohammed Elwardi Fadeli (elwardifadeli@gmail.com)

\*---------------------------------------------------------------------------*/

#ifndef testMacros_H
#define testMacros_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include <mpi.h>
#include "Pstream.H"
#include "memberStealer.H"

// Time-related macros

// Specialize the member stealer for case and processorCase
#define prepareTimePaths() \
   SPECIALIZE_MEMBER_STEALER(case_, fileName, TimePaths); \
   SPECIALIZE_MEMBER_STEALER(processorCase_, bool, TimePaths)

// Reset time case name and processorCase
#define modifyTimePaths(timeName, newProcessorCase, newCaseName) \
    GET_REF(processorCase_, TimePaths, timeName) = newProcessorCase; \
    GET_REF(case_, TimePaths, timeName) = newCaseName

// Mesh-related macros

// Reset mesh pointer Sample dynamic Mesh object
#define resetMeshPointer(timeObj, meshPtr, meshClass, region)      \
meshPtr.reset \
( \
    new meshClass \
    ( \
        IOobject \
        ( \
            region, \
            timeObj.timeName(), \
            timeObj, \
            Foam::IOobject::MUST_READ \
        ) \
    ) \
)

// Field-related macros

// Create a field from an IOobject, mesh, and initial value
#define createField(fieldName, fieldClass, time, mesh, dimSet, zero) \
fieldClass fieldName \
( \
    IOobject \
    ( \
        #fieldName, \
        time.timeName(), \
        mesh, \
        IOobject::NO_READ, \
        IOobject::AUTO_WRITE \
    ), \
    mesh, \
    dimensionedScalar(#fieldName, dimSet, zero) \
)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
