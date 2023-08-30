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

\*---------------------------------------------------------------------------*/

#ifndef MultiSpeciesTransportModel_H
#define MultiSpeciesTransportModel_H

#include "multiSpeciesTransportModel.H"
#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseModel;
class phaseSystem;

/*---------------------------------------------------------------------------*\
                           Class MultiSpeciesTransportModel Declaration
\*---------------------------------------------------------------------------*/
template<class Thermo>
class MultiSpeciesTransportModel
:
    public multiSpeciesTransportModel
{

protected:

        //- Thermo
        const Thermo& thermo_;

public:

    // Constructors

        MultiSpeciesTransportModel
        (
            const phaseModel& phase
        );

    //- Destructor
        virtual ~MultiSpeciesTransportModel();

    // Member Functions

        //- Access function to incompressible transport model
    virtual const Thermo& thermo() const
    {
        return thermo_;
    }

        //- Effective thermal turbulent diffusivity for temperature
    //  of mixture [W/m/K]
    virtual tmp<volScalarField> kappaEff() const
    {
        return thermo_.kappa();
    }

    // //- Effective thermal turbulent diffusivity for temperature
    // //  of mixture for patch [W/m/K]
    // virtual tmp<scalarField> kappaEff(const label patchi) const
    // {
    //     return thermo_.kappa().boundaryField()[patchi];
    // }

    //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
    virtual tmp<volScalarField> alphaEff() const
    {
        return thermo_.alphahe();
    }

    //- Effective thermal turbulent diffusivity of mixture
    //  for patch [kg/m/s]
    virtual tmp<scalarField> alphaEff(const label patchi) const
    {
        return thermo_.alphahe(patchi);
    }



};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MultiSpeciesTransportModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
