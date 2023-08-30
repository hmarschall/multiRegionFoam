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
    Foam::electricConductive.C

SourceFiles
    electricConductive.C

\*---------------------------------------------------------------------------*/

#ifndef electricConductive_H
#define electricConductive_H

#include "regionType.H"
#include "dimensionedScalar.H"
#include "multiRegionSystem.H"
#include "volFields.H"
#include "IOReferencer.H"
#include "HashPtrTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class sigmaModelList;
class dissolvedModel;

namespace regionTypes
{

/*---------------------------------------------------------------------------*\
                           Class electricConductive Declaration
\*---------------------------------------------------------------------------*/

class electricConductive
:
    public regionType
{
        //- name
        word regionName_;

        //- transport properties dictionary 
        IOdictionary electricProperties_;

        //- Current density (A/m2)
        volVectorField i_;

        //- Potential (V)
        volScalarField phi_;

        //- Source term (A/m3)
        volScalarField j_;

        //- Electric conductivity
        volScalarField sigmaField_;

        //- Temperature
        volScalarField T_;

        //- Sigma model (electric conductivity)
        autoPtr<sigmaModelList> sigma_;

        //- Relaxation
        scalar relax_;

        //- Active side?
        //- Reference side: voltage = 0.0
        //- Active side: voltage: fixed or calculated
        Switch active_;

        //- Dissolved water
        //- In PEM fuel cell/electrolyzer, water can dissolve from one side
        //- to the other. This process is usually affected by the proton transfer.
        //- Refer to the region of CCM (carbon coated membrane)
        Switch dissolveOnOff_;
        autoPtr<dissolvedModel> dissolved_;

        //- Patch name
        //- The reference boundary for the prescribed current density or voltage
        word patchName_;

        //- Operation mode: (currentCollector_ true)
        Switch galvanostatic_;

        //- Mean current density
        //- Make sense when galvanostatic_ is true
        scalar ibar_;

        //- Voltage
        //- Make semse when galvanostatic_ is false
        scalar voltage_;

        IOdictionary cellProperties_;

protected:



public:

    TypeName("electricConductive");

    // Constructors

        electricConductive
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~electricConductive();


    // Member Functions

        // Access


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

        virtual scalar getMinDeltaT();

        virtual void solveRegion();

        virtual void prePredictor();

        virtual void momentumPredictor();

        virtual void pressureCorrector();

        virtual void meshMotionCorrector();

        virtual void setCoupledEqns();

        virtual void postSolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
