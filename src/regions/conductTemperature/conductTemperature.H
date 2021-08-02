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
public:

    typedef
        HashTable
        <
            Foam::tmp<Foam::fvScalarMatrix>, word, word::hash
        >
        coupledScalarEqns;

private:

        //- mesh
        const fvMesh& mesh_;

        //- name
        word regionName_;

        //- transport properties dictionary 
        IOdictionary transportProperties_;

        //- thermal conductivity
        dimensionedScalar k_;

        //- specific heat
        dimensionedScalar cv_;

        //- density
        dimensionedScalar rho_;

        //- temperature
        volScalarField T_;

        //- table of coupled equations
        coupledScalarEqns coupledScalarEqnsTable_;

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

        // coupledScalarEqns coupledScalarEqnsTable() const;

        virtual void correct();

        virtual void setRDeltaT();

        virtual void solveRegion();

        virtual void solveCoupledPartitioned();

        // virtual tmp<fvScalarMatrix> coupledFvScalarMatrix() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
