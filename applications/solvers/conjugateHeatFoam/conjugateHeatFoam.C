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

Application
    conjugateHeatPhaseFieldFoam

Description
    Multi-region conjugate heat transfer solver for thermal energy transport 
    between two separate regions (fluid and solid) supporting Dirichlet-Neumann 
    monolithic and partitioned coupling methods.
    
    Additionally, the solver permits a thermal contact resistance between 
    regions via a "resistive" interface boundary condition.

Author
    Holger Marschall (holger.marschall@tu-darmstadt.de). All rights reserved.
    Brent A. Craven (craven@psu.edu). All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "coupledFvMatrices.H"

#include "subCycle.H"

#include "regionCoupleTemperatureFvPatchScalarField.H"
#include "regionCoupleHeatFluxFvPatchScalarField.H"
#include "regionCouplePolyPatch.H"

//- Supported coupling methods
enum couplingMethods
{
    partitioned,
    monolithic
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "readCoupledSolutionControls.H"
    #include "createFields.H"

    DynamicList<label> fluidCoupledTemperaturePatchIDs(T.boundaryField().size());
    DynamicList<label> fluidCoupledHeatFluxPatchIDs(T.boundaryField().size());

    DynamicList<label> solidCoupledTemperaturePatchIDs(Ts.boundaryField().size());
    DynamicList<label> solidCoupledHeatFluxPatchIDs(Ts.boundaryField().size());

    if (couplingMethod == partitioned)
    {
        #include "checkCoupledPatches.H"
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        // Detach patches
        if (couplingMethod == monolithic)
        {
            #include "detachPatches.H"
        }

        runTime++;

        Info<< nl << "Time = " << runTime.timeName() << nl << endl;

        // Coupled patches
        if (couplingMethod == monolithic)
        {
            #include "attachPatches.H"
        }

        #include "TEqns.H"

        if (runTime.outputTime())
        {
            runTime.write();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << "s" << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
