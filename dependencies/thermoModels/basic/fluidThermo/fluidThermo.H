/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fluidThermo

Description
    Basic thermodynamic properties based on density

SourceFiles
    fluidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef fluidThermo_H
#define fluidThermo_H

#include "basicThermo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class fluidThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidThermo
:
    public basicThermo
{

public:

    //- Runtime type information
    TypeName("fluidThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // //- Declare run-time constructor selection table for dictionary based
    // declareRunTimeSelectionTable
    // (
    //     autoPtr,
    //     fluidThermo,
    //     fvMeshDictPhase,
    //     (const fvMesh& mesh, const word& phaseName, const word& dictName),
    //     (mesh, phaseName, dictName)
    // );

    // Constructors

        //- Construct from mesh and phase name
        fluidThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        // //- Construct from mesh and phase name
        // fluidThermo
        // (
        //     const fvMesh&,
        //     const dictionary&,
        //     const word& phaseName
        // );


        // //- Construct from mesh,dictionary,phase name with a single temperature
        // fluidThermo
        // (
        //     const fvMesh&,
        //     const word& phaseName,
        //     const word& dictionaryName
        // );


    //- Selector
    static autoPtr<fluidThermo> New
    (
        const fvMesh&,
        const word& phaseName=word::null
    );


    // //- Selector
    // static autoPtr<fluidThermo> New
    // (
    //     const fvMesh&,
    //     const word& phaseName,
    //     const word& dictName
    // );


    //- Destructor
    virtual ~fluidThermo();


    // Member functions

        // Access to thermodynamic state variables

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            //  Limit thermo rho between rhoMin and rhoMax
            virtual void correctRho
            (
                const volScalarField& deltaRho,
                const dimensionedScalar& rhoMin,
                const dimensionedScalar& rhoMax
            ) = 0;


            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            virtual void correctRho
            (
                const volScalarField& deltaRho
            ) = 0;

            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const = 0;


        // Access to transport state variables

            // TODO: See basicThermo
            //- Dynamic viscosity of mixture [kg/m/s]
            virtual tmp<volScalarField> mu() const = 0;

            //- Dynamic viscosity of mixture for patch [kg/m/s]
            virtual tmp<scalarField> mu(const label patchi) const = 0;

            //- Kinematic viscosity of mixture [m^2/s]
            virtual tmp<volScalarField> nu() const;

            //- Kinematic viscosity of mixture for patch [m^2/s]
            virtual tmp<scalarField> nu(const label patchi) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
