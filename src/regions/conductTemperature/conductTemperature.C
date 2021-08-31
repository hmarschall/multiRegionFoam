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

#include "fvCFD.H"
#include "conductTemperature.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

#include "simpleControl.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(conductTemperature, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        conductTemperature,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::conductTemperature::conductTemperature
(
    const fvMesh& mesh,
    const word& regionName
)
:
    regionType(mesh, regionName),

    mesh_(mesh),
    regionName_(regionName),

    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            this->time().constant(),
            *this,
//            this->time().timeName(),
//            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    k_(transportProperties_.lookup("k")),
    cv_(transportProperties_.lookup("cv")),
    rho_(transportProperties_.lookup("rho")),
    T_
    (
        IOobject
        (
            "T",
            this->time().timeName(),
            *this,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        *this
    )
{
//    eqns_.insert
//    (
//        T_.name() + "Eqn",
//        new fvScalarMatrix
//        (
//            T_,
//            dimensionSet(1, -1, -2, 1, 0, 0, 0)
//        )
//    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::conductTemperature::~conductTemperature()
{
//    delete TEqnPtr_;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::conductTemperature::correct()
{
    // do nothing, add as required
}


void Foam::regionTypes::conductTemperature::setRDeltaT()
{
    // do nothing, add as required
}


void Foam::regionTypes::conductTemperature::setCoupledEqns()
{
    fvScalarMatrix TEqn =
    (
        fvm::ddt(rho_*cv_, T_)
     ==
        fvm::laplacian(k_, T_, "laplacian(k,T)")
    );

    fvScalarMatrices.set
    (
        T_.name() + this->name() + "Eqn",
        new fvScalarMatrix(TEqn)
    );

//    TEqn_.relax();
//    TEqn_.solve();



//    TEqnPtr_ =
//        new IOReferencer<fvScalarMatrix>
//        (
//            IOobject
//            (
//                T_.name() + "Eqn",
//                this->time().timeName(),
//                *this,
//                IOobject::NO_READ,  /*must be NO_READ*/
//                IOobject::NO_WRITE  /*must be NO_WRITE*/
//            ),
//            &TEqn_
////            solidInterface::New(solIntType, sigma.mesh(), *this).ptr()
//        );


//    regCoupledEqn(TEqn_);

//    TEqn.relax();
//    TEqn.solve();

//    TEqnPtr_ = new fvScalarMatrix
//    (
//        fvm::ddt(rho_*cv_, T_)
//     ==
//        fvm::laplacian(k_, T_, "laplacian(k,T)")
//    );

//    TEqnPtr_->solve();

//    regCoupledEqn(*TEqnPtr_);
}

void Foam::regionTypes::conductTemperature::solveRegion()
{
    // do nothing, add as required

//    Info << nl << "Solving for temperature in " << regionName_ << endl;
//    simpleControl simpleControlRegion(*this);

//    while (simpleControlRegion.correctNonOrthogonal())
//    {
//        TEqn->relax();
//        TEqn->solve();
//    }
}

//template<class Type>
//void Foam::regionTypes::conductTemperature::regCoupledEqn
//(
//    fvMatrix<Type>& fvm
////    const DimensionedField<Type, volMesh>& psi
//)
//{
////    const fvMesh mesh = psi.mesh();

//    // checkin to object registry if not already present
////    if
////    (
////       !(
////            this->thisDb().foundObject<IOReferencer<fvMatrix<Type> > >
////            (
////                fvm.psi().name() + "Eqn"
////            )
////        )
////    )
//    {
//        IOReferencer<fvMatrix<Type> >* fvmRef =
//        new IOReferencer<fvMatrix<Type> >
//        (
//            IOobject
//            (
//                fvm.psi().name() + "Eqn",
//                this->time().timeName(),
//                *this,
//                IOobject::NO_READ,  /*must be NO_READ*/
//                IOobject::NO_WRITE  /*must be NO_WRITE*/
//            ),
//            &fvm
////            &(const_cast<fvMatrix<Type>& >(fvm))
//        );

//        Info<< "Registered " << fvmRef->name() 
//            << " from region " << this->name()
//            << " to object registry" << nl << endl;
//    }
//}


// ************************************************************************* //
