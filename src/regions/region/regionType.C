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

\*---------------------------------------------------------------------------*/

#include "regionType.H"
#include "multiRegionSystem.H"
#include "IOReferencer.H"

namespace Foam
{
    defineTypeNameAndDebug(regionType, 0);
    defineRunTimeSelectionTable(regionType, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::regionType::regionType
(
    const fvMesh& mesh,
    const word& regionName
)
:
    fvMesh
    (
        IOobject
        (
            regionName,
            mesh.time().timeName(),
            mesh.time(),
            IOobject::MUST_READ
        )
    ),

    mesh_(mesh),

    dict_
    (
        IOobject
        (
            "multiRegionProperties",
            mesh_.time().constant(),
            mesh_.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//-HM: see compressible/steadyUniversalMRFFoam/createThermo.H  43   
// psisPtr = const_cast<volScalarField*>(&(thermoPtr->psi()));

//template<class Type>
//void Foam::regionType::regCoupledEqn
//(
//    const fvMatrix<Type>& fvm,
//    const fvMesh& mesh
//)
//{
//    //const fvMesh& mesh = fvm.psi().mesh(); //! not working (private in this context)
//
//    // checkin to object registry if not already present
//    if
//    (
//       !(
//            mesh.thisDb().foundObject<IOReferencer<fvMatrix<Type> > >
//            (
//                fvm.psi().name() + "Eqn"
//            )
//        )
//    )
//    {
//        IOReferencer<fvMatrix<Type> >* fvmRef =
//        new IOReferencer<fvMatrix<Type> >
//        (
//            IOobject
//            (
//                fvm.psi().name() + "Eqn",
//                mesh.time().timeName(),
//                mesh,
//                IOobject::NO_READ,  /*must be NO_READ*/
//                IOobject::NO_WRITE  /*must be NO_WRITE*/
//            ),
//            &(const_cast<fvMatrix<Type>& >(fvm))
//        );

//        Info<< "Registered " << fvmRef->name() 
//            << " from region " << mesh.name()
//            << " to object registry" << nl << endl;
//    }
//}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionType::~regionType()
{}

// ************************************************************************* //
