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
    Foam::conductTemperature.C

SourceFiles
    conductTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef conductTemperature_H
#define conductTemperature_H

#include "regionType.H"
#include "dimensionedScalar.H"
#include "multiRegionSystem.H"
#include "volFields.H"
#include "IOReferencer.H"
#include "HashPtrTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{

/*---------------------------------------------------------------------------*\
                           Class conductTemperature Declaration
\*---------------------------------------------------------------------------*/

class conductTemperature
:
    public regionType
{
        //- mesh
        const fvMesh& mesh_;

        //- name
        word regionName_;

        //- transport properties dictionary 
        IOdictionary transportProperties_;

        //- specific heat
        // dimensionedScalar cv_;
        volScalarField cv_;

        //- density
        // dimensionedScalar rho_;
        volScalarField rho_;

        //- thermal diffusivity
        autoPtr<volScalarField> alpha_;

        //- thermal conductivity
        autoPtr<volScalarField> k_;

protected:

        //- temperature
//        volScalarField T_;

        autoPtr<volScalarField> T_;

public:

    TypeName("conductTemperature");

    // Constructors

        conductTemperature
        (
            const fvMesh& mesh,
            const word& regionName
        );

    //- Destructor
    virtual ~conductTemperature();


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
                this->thisDb().lookupObject<T>(name)
            );
        }

        template<class T>
        T& getObject(const word& name)
        {
            return
            (
                const_cast<const T&>(this->thisDb().lookupObject<T>(name))
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
        
        virtual void updateFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
