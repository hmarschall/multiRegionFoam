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
#include "transportTemperature.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

#include "simpleControl.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(transportTemperature, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        transportTemperature,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::transportTemperature::transportTemperature
(
    const fvMesh& mesh,
    const word& regionName
)
:
    regionType(mesh, regionName),

    mesh_(mesh),
    regionName_(regionName),

    U_
    (
        IOobject
        (
            "U",
            this->time().timeName(),
            *this,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        *this
    ),
    phi_
    (
		IOobject
		(
			"phi",
            this->time().timeName(),
            *this,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		linearInterpolate(U_) & (*this).Sf()    
    ),

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
    cp_(transportProperties_.lookup("cp")),
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
    ),
    TEqnPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::transportTemperature::~transportTemperature()
{
    delete TEqnPtr_;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::transportTemperature::correct()
{
    // do nothing, add as required
}


void Foam::regionTypes::transportTemperature::setRDeltaT()
{
    // do nothing, add as required
}


void Foam::regionTypes::transportTemperature::solveCoupledPartitioned()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportTemperature::solveRegion()
{
    Info << nl << "Solving for temperature in " << regionName_ << endl;
    simpleControl simpleControlRegion(*this);

    dimensionedScalar alpha = k_/(rho_*cp_);

    while (simpleControlRegion.correctNonOrthogonal())
    {
        TEqnPtr_ = new fvScalarMatrix
        (
            fvm::ddt(T_)
          + fvm::div(phi_, T_)
         ==
            fvm::laplacian(alpha, T_)
        );

//        regCoupledEqn(TEqn, T_.mesh());
        regCoupledEqn(*TEqnPtr_);

        TEqnPtr_->relax();
        TEqnPtr_->solve();
    }


//    dimensionedScalar alpha = k_/(rho_*cp_);

//    fvScalarMatrix* TEqn = new fvScalarMatrix
//    (
//        fvm::ddt(T_)
//      + fvm::div(phi_, T_)
//      ==
//        fvm::laplacian(alpha, T_)
//    );

//    IOReferencer<fvScalarMatrix> TEqnRef
//    (
//        IOobject
//        (
//            TEqn->psi().name() + "Eqn",
//            this->time().timeName(),
//            *this,
//            IOobject::NO_READ,  /*must be NO_READ*/
//            IOobject::NO_WRITE  /*must be NO_WRITE*/
//        ),
//        TEqn
//    );

//    TEqn->relax();
//    TEqn->solve();


//    const objectRegistry& db = *this; //this->thisDb();

//    Info << nl <<"Objects registered to "
//         << this->name() << " :" << nl << endl;

//    forAllConstIter(HashTable<regIOobject*>, db, iter)
//    {
//        Info << " name : " << iter()->name() << nl
//             << " type : " << iter()->type() << nl
//             << endl;
//    }

}

template<class Type>
void Foam::regionTypes::transportTemperature::regCoupledEqn
(
    const fvMatrix<Type>& fvm
//    const DimensionedField<Type, volMesh>& psi
)
{
//    const fvMesh mesh = psi.mesh();

    // checkin to object registry if not already present
    if
    (
       !(
            this->thisDb().foundObject<IOReferencer<fvMatrix<Type> > >
            (
                fvm.psi().name() + "Eqn"
            )
        )
    )
    {
        IOReferencer<fvMatrix<Type> >* fvmRef =
        new IOReferencer<fvMatrix<Type> >
        (
            IOobject
            (
                fvm.psi().name() + "Eqn",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,  /*must be NO_READ*/
                IOobject::NO_WRITE  /*must be NO_WRITE*/
            ),
            &(const_cast<fvMatrix<Type>& >(fvm))
        );

        Info<< "Registered " << fvmRef->name() 
            << " from region " << this->name()
            << " to object registry" << nl << endl;
    }
}

// tmp<fvScalarMatrix> 
// Foam::regionTypes::transportTemperature::coupledFvScalarMatrix() const
// {
//     dimensionedScalar alpha = k_/(rho_*cp_);

//     return tmp<fvScalarMatrix>
//     (
//         new fvScalarMatrix
//         (
//             fvm::ddt(T_)
//         + fvm::div(phi_, T_)
//         ==
//             fvm::laplacian(alpha, T_)
//         )
//     );
// }

// ************************************************************************* //
