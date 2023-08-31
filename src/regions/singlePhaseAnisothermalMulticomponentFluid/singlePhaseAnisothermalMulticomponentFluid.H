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
    Foam::singlePhaseAnisothermalMulticomponentFluid.C

SourceFiles
    singlePhaseAnisothermalMulticomponentFluid.C

Author
    Steffen Hess(s.hess@fz-juelich.de, Affiliation A) 
    All rights reserved.

Contact
    Steffen Hess(s.hess@fz-juelich.de) 

Affiliations
    Affiliation A)
    Forschungszentrum Jülich IEK-14

Acknowledgement
    Funded by
    AIDAS – AI, Data Analytics and Scalable Simulation

Description
    
    This file is part of the multiRegionFoam library.

    You may refer to this software as :
    https://doi.org/10.48550/arXiv.2306.01924

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#ifndef singlePhaseAnisothermalMulticomponentFluid_H
#define singlePhaseAnisothermalMulticomponentFluid_H

#include "regionType.H"
#include "dimensionedScalar.H"
#include "multiRegionSystem.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOReferencer.H"
#include "HashPtrTable.H"
#include "zeroGradientFvPatchField.H"
#include "pimpleControl.H"
// #include "singlePhaseTransportModel.H"
// #include "RASModel.H"
#include "initContinuityErrs.H"

#include "rhoReactionThermo.H"
#include "multivariateScheme.H"
#include "phaseSystem.H"
#include "phaseModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{

/*---------------------------------------------------------------------------*\
                           Class singlePhaseAnisothermalMulticomponentFluid Declaration
\*---------------------------------------------------------------------------*/

class singlePhaseAnisothermalMulticomponentFluid
:
    public regionType
{

        //- name
        word regionName_;

        //- phase type
        autoPtr<phaseSystem> phases_;

        //- Cached temperature matrix
        tmp<fvScalarMatrix> TEqn;


protected:

// TODO: I do not see this, where??
// Why we are using autoPtr for the fields
// Weil bei der Instantiierung der Fields bei gekoppelter RB bereits das benachbarte Feld benötigt wird,
// dieses aber ggf. noch nicht da ist, weil die region später instantiiert wird. 
// Jetzt kann man den autoPtr im Konstruktor anlegen und später auf ein Feld zeigen lassen.

        //- temperature
//        volScalarField T_;

        // autoPtr<volScalarField> T_;

public:

    TypeName("singlePhaseAnisothermalMulticomponentFluid");

    // Constructors

        singlePhaseAnisothermalMulticomponentFluid
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~singlePhaseAnisothermalMulticomponentFluid();


    // Member Functions

        // // Access
        // const volScalarField& T() const
        // {
        //     return T_();
        // }

        // volScalarField& T()
        // {
        //     return T_();
        // }

        word regionTypeName() const
        {
            return singlePhaseAnisothermalMulticomponentFluid::typeName;
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
