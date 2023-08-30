/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
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
    Foam::phaseSystem

Description
    Class to represent a system of phases and model interfacial transfers
    between them.

SourceFiles
    phaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef phaseSystem_H
#define phaseSystem_H

#include "IOdictionary.H"

#include "phaseModel.H"
// #include "phaseInterface.H"
// #include "phaseInterfaceKey.H"
#include "HashPtrTable.H"
// #include "PtrListDictionary.H"
#include "PtrList.H"
#include "hashedWordList.H"
#include "pimpleControl.H"
// #include "fvOptions.H"

// #include "IOMRFZoneList.H"
// #include "fvModels.H"
// #include "fvConstraints.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"

#include "porosityModelList.H"

#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// class surfaceTensionModel;
// class pressureReference;
// class nonOrthogonalSolutionControl;

class pimpleControl;

/*---------------------------------------------------------------------------*\
                         Class phaseSystem Declaration
\*---------------------------------------------------------------------------*/

class phaseSystem
:
    public IOdictionary
{
public:

        // typedef PtrListDictionary<phaseModel> phaseModelList;
        typedef PtrList<phaseModel> phaseModelList;


protected:
    // Protected data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Phase models
        phaseModelList phaseModels_;

        //- Porous zone list
        autoPtr<porosityModelList> porosityModels_;

        //- Pressure p_rgh
        volScalarField p_rgh_;

        //- pRefCell
        label pRefCell_;

        //- pRefValue
        scalar pRefValue_;

        //- Pimple control
        pimpleControl pimple_;

public:

    //- Runtime type information
    TypeName("phaseSystem");

    //- Default name of the phase properties dictionary
    static const word propertiesName;


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseSystem,
            dictionary,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );


    // Constructors

        //- Construct from fvMesh
        phaseSystem(const fvMesh& mesh);


    // Selectors

        static autoPtr<phaseSystem> New
        (
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~phaseSystem();


    // Member Functions

        void reset(const IOdictionary& dict);

        // Access

            //- Return the mesh
            inline const fvMesh& mesh() const;

            //- Return the phase models
            inline const phaseModelList& phases() const;

            //- Access the phase models
            inline phaseModelList& phases();

            inline const porosityModelList& porosityModels() const;

            inline porosityModelList& porosityModels();

            // //- Access the fvOptions
            // inline fv::options& fvOptions() const;

        // Evolution

            //- Solve for the phase fractions
            virtual void solve() = 0;

            //- Correct the fluid properties other than those listed below
            virtual void correct();

            //- Correct the thermodynamics
            virtual void correctThermo();

            //- Correct the energy transport e.g. alphat
            virtual void correctEnergyTransport();

        // Surface coupling

            //- Return the TEqn from the phases
            virtual tmp<fvScalarMatrix> TEqn() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseSystemI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
