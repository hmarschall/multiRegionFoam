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
    Foam::transportTemperature.C

SourceFiles
    transportTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef transportTemperature_H
#define transportTemperature_H

#include "regionType.H"
#include "dimensionedScalar.H"
#include "multiRegionSystem.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOReferencer.H"
#include "HashPtrTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{

/*---------------------------------------------------------------------------*\
                           Class transportTemperature Declaration
\*---------------------------------------------------------------------------*/

class transportTemperature
:
    public regionType
{
        //- mesh
        const fvMesh& mesh_;

        //- name
        word regionName_;

        //- transport
        volVectorField U_;

        surfaceScalarField phi_;

        //- transport properties dictionary 
        IOdictionary transportProperties_;

        //- thermal conductivity
        dimensionedScalar k_;

        //- specific heat
        dimensionedScalar cp_;

        //- density
        dimensionedScalar rho_;

protected:

        //- temperature
        volScalarField T_;

        //- equation-equation coupling
//        template<class Type>
//        struct eqnTable
//        {
//            HashPtrTable<Type> eqns_;
//        };

//        eqnTable<scalar> eqns_;

//        HashPtrTable<fvScalarMatrix> eqns_;

//        IOReferencer<fvScalarMatrix>* TEqnPtr_;

public:

    TypeName("transportTemperature");

    // Constructors

        transportTemperature
        (
            const fvMesh& mesh,
            const word& regionName
        );

    //- Destructor
    virtual ~transportTemperature();


    // Member Functions

        // Access
        const volScalarField& T() const
        {
            return T_;
        }

        volScalarField& T()
        {
            return T_;
        }

        virtual bool active() const
        {
            return true;
        }

        virtual void correct();

        virtual void setRDeltaT();

        virtual void solveRegion();

        virtual void setCoupledEqns();

        // Get equation by name
//        fvScalarMatrix& getCoupledEqn
//        (
//            word name
//        )
//        {
//            if (!fvScalarMatrices.table.found(name))
//            {
//                FatalErrorIn("transportTemperature::getCoupledEqn")
//                    << "Equation " << name << " does not exist."
//                    << abort(FatalError);
//            }

//            return *fvScalarMatrices.table[name];
//        }

//        template<class Type>
//        fvMatrix<Type>& getCoupledEqn
//        (
//            word name
//        )
//        {
//            if (Type::componentName == "scalar")
//            {
//                if (!fvScalarMatrices.table.found(name))
//                {
//                    FatalErrorIn("transportTemperature::getCoupledEqn")
//                        << "Equation " << name << " does not exist."
//                        << abort(FatalError);
//                }

//                return *fvScalarMatrices.table[name];
//            }
//            else if (Type::componentName == "vector")
//            {
//                if (!fvVectorMatrices.table.found(name))
//                {
//                    FatalErrorIn("transportTemperature::getCoupledEqn")
//                        << "Equation " << name << " does not exist."
//                        << abort(FatalError);
//                }

//                return *fvVectorMatrices.table[name];
//            }
//            else if (Type::componentName == "symmTensor")
//            {
//                if (!fvSymmTensorMatrices.table.found(name))
//                {
//                    FatalErrorIn("transportTemperature::getCoupledEqn")
//                        << "Equation " << name << " does not exist."
//                        << abort(FatalError);
//                }

//                return *fvSymmTensorMatrices.table[name];
//            }            
//            else if (Type::componentName == "tensor")
//            {
//                if (!fvTensorMatrices.table.found(name))
//                {
//                    FatalErrorIn("transportTemperature::getCoupledEqn")
//                        << "Equation " << name << " does not exist."
//                        << abort(FatalError);
//                }

//                return *fvTensorMatrices.table[name];
//            }
//            else
//            {
//                FatalErrorIn("transportTemperature::getCoupledEqn")
//                    << "Equation type not supported"
//                    << abort(FatalError);
//            }

//            return (nullptr);
//        }

//        template<class Type>
//        void regCoupledEqn
//        (
//            fvMatrix<Type>& fvm
//        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
