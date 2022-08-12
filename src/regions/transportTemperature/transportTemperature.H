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
        //- name
        word regionName_;

        //- transport
        //volVectorField U_;

        //surfaceScalarField phi_;

        //- transport properties dictionary 
        IOdictionary transportProperties_;

        //- thermal conductivity
        dimensionedScalar k_;

        //- specific heat
        dimensionedScalar cp_;

        //- density
        dimensionedScalar rho_;

        //- transport
        autoPtr<volVectorField> U_;

        //- thermal diffusivity
        autoPtr<volScalarField> alpha_;

        //- flux
        autoPtr<surfaceScalarField> phi_;

protected:

        //- temperature
        autoPtr<volScalarField> T_;

public:

    TypeName("transportTemperature");

    // Constructors

        transportTemperature
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~transportTemperature();


    // Member Functions

        // Access
        const volScalarField& T() const
        {
            return T_();
        }

        volScalarField& T()
        {
            return T_();
        }

        // Access by name
        template<class T>
        const T& getObject(const word& name)
        {
            return
            (
                mesh().thisDb().lookupObject<T>(name)
            );
        }

        template<class T>
        T& getObject(const word& name)
        {
            return
            (
                const_cast<const T&>(mesh().thisDb().lookupObject<T>(name))
            );
        }

        virtual bool active() const
        {
            return true;
        }

        virtual void correct();

        virtual void setRDeltaT();

        virtual void solveRegion();

        virtual void setCoupledEqns();
        
        virtual void postSolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
