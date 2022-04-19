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
    Foam::singlePhaseAnisothermalMulticomponentFluid.C

SourceFiles
    singlePhaseAnisothermalMulticomponentFluid.C

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

#include "hReactionThermo.H"
#include "multivariateScheme.H"
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
        //- global mesh
        const fvMesh& mesh_;

        //- name
        word regionName_;

        //- phase dictionary
        const IOdictionary dict_;

        // //- thermoModel
        autoPtr<hReactionThermo> thermo_;

        //- phaseModel, including the thermo model
        autoPtr<phaseModel> electrochemicalPhaseModel_;

        //- Pimple control
        pimpleControl pimple_;

        //- Pressure p_rgh
        volScalarField p_rgh_;

protected:

        //- temperature
//        volScalarField T_;

        // autoPtr<volScalarField> T_;

public:

    TypeName("singlePhaseAnisothermalMulticomponentFluid");

    // Constructors

        singlePhaseAnisothermalMulticomponentFluid
        (
            const fvMesh& mesh,
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
