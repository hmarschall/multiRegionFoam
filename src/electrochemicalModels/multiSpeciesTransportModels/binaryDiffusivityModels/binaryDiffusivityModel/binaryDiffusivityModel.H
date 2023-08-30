/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

Namespace
    Foam::binaryDiffusivityModels

Description
    Namespace for binary diffusivity models.

Class
    Foam::binaryDiffusivityModel

Description
    An abstract base class for binary diffusivity models.

SourceFiles
    binaryDiffusivityModel.C

\*---------------------------------------------------------------------------*/

#ifndef binaryDiffusivityModel_H
#define binaryDiffusivityModel_H

#include "IOdictionary.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"
#include "volFieldsFwd.H"
#include "phaseSystem.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class binaryDiffusivityModel Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class binaryDiffusivityModel
{

protected:

    // Protected data
    
        //- Transport properties IOdictionary
        const IOdictionary& dict_;

        //- Species names
        const word name1_;

        const word name2_;

        const ThermoType& thermo_;

        //- phaseModel
        const phaseModel& phase_;

        //- phaseSystem
        const phaseSystem& phaseSys_;

        //- Porous zones
        const porosityModelList& pZones_ ;

    // Private Member Functions

        //- Disallow copy construct
        binaryDiffusivityModel(const binaryDiffusivityModel&);

        //- Disallow default bitwise assignment
        void operator=(const binaryDiffusivityModel&);


public:

    // //- Runtime type information
    // TypeName("binaryDiffusivityModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            binaryDiffusivityModel,
            dictionary,
            (
                const ThermoType& thermo,
                const word name1,
                const word name2,
                const IOdictionary& dict,
                const phaseModel& phase
            ),
            (thermo, name1, name2, dict, phase)
        );


    // Selectors

        //- Return a reference to the selected diffusivity model
        static autoPtr<binaryDiffusivityModel> New
        (
            const ThermoType& thermo,
            const word name1,
            const word name2,
            const IOdictionary& dict,
            const phaseModel& phase
        );


    // Constructors

        //- Construct from components
        binaryDiffusivityModel
        (
            const ThermoType& thermo,
            const word name1,
            const word name2,
            const IOdictionary& dict,
            const phaseModel& phase
        );


    // Destructor

        virtual ~binaryDiffusivityModel()
        {
            Info << "Calling Destructor of binaryDiffusivityModel " << endl;
        }


    // Member Functions

         //- Binary diffusion coefficient [m^2/s]
         virtual tmp<volScalarField> D() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// #ifdef NoRepository
//     #include "newBinaryDiffusivityModel.C"
// #endif

#ifdef NoRepository
    #include "binaryDiffusivityModel.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
