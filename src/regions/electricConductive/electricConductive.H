/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Fork:     foam-extend
    \\  /    A nd           | Version:  4.1                                 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::electricConductive.C

SourceFiles
    electricConductive.C

Author
    Steffen Hess(s.hess@fz-juelich.de, Affiliation A) 
    Constantin Habes (constantin.habes@tu-darmstadt.de, Affiliation B)
    Holger Marschall (holger.marschall@tu-darmstadt.de, Affiliation B)

Contact
    Steffen Hess(s.hess@fz-juelich.de) 
    Holger Marschall (holger.marschall@tu-darmstadt.de)
    main developer and principal investigator
    TU Darmstadt

Affiliations
    Affiliation A)
    Forschungszentrum Jülich IEK-14
    Affiliation B)
    Computational Multiphase Flow
    Department of Mathematics
    Technical University of Darmstadt, Germany

Acknowledgement
    Funded by
    AIDAS – AI, Data Analytics and Scalable Simulation and
    Hessian Ministry of Higher Education, Research, Science and the Arts,
    National High Performance Computing Center for Computational Engineering
    Science (NHR4CES) and 

Description
    
    This file is part of the multiRegionFoam library.

    You may refer to this software as :
    https://doi.org/10.48550/arXiv.2306.01924

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

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
