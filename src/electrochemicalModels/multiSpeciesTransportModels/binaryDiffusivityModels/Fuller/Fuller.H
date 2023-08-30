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

Class
    Foam::binaryDiffusivityModels::Fuller

Description
    Fuller binary diffusivity model.
    
SourceFiles
    Fuller.C

\*---------------------------------------------------------------------------*/

#ifndef Fuller_H
#define Fuller_H

#include "binaryDiffusivityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "fsgMolecularWeights.H"
#include "fsgDiffusionVolumes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// namespace binaryDiffusivityModels
// {

/*---------------------------------------------------------------------------*\
                           Class Fuller Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class Fuller
:
    public binaryDiffusivityModel<ThermoType>
{
    // Private data

        const dictionary binaryDiffDict_;

        scalar W12;     // [g/mol]

        scalar V12;     // [cm^3]
      
        //- Conversion factor to convert the binary diffusivity coefficient
        //  in standard units (m^2/s)
        scalar pAtm;

        // Units correction factor
        dimensionedScalar unitsCorr_;      


public:

    //- Runtime type information
    TypeName("Fuller");


    // Constructors

        //- construct from components
        Fuller
        (
            const ThermoType& thermo,
            const word name1,
            const word name2,
            const IOdictionary& dict,
            const phaseModel& phase
        );


    // Destructor

        ~Fuller()
        {}


    // Member Functions

         //- Binary diffusion coefficient [m^2/s]
         virtual tmp<volScalarField> D() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// } // End namespace binaryDiffusivityModels
} // End namespace Foam

#ifdef NoRepository
    #include "Fuller.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
