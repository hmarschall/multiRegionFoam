/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    updateMesh

Description


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "dynamicFvMesh.H"
#include "dynamicRefineFvMesh.H"

#include "ReadFields.H"
#include "IOobjectList.H"

#include <iostream>
#include <fstream>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    argList::validArgs.clear();
#   include "addTimeOptions.H"
    argList::validOptions.insert("overwrite", "");

#   include "setRootCase.H"
#   include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createDynamicFvMesh.H"

    bool overwrite = args.optionFound("overwrite");

#   include "createFields.H"

//    for (label i = 0; i < refineInterval; i++)
//    {
//        if (!overwrite)
//        {
//            runTime++;
//        }

//        mesh.update();

//        Info << "Writing modified mesh to time " << runTime.timeName() << endl;
//    }

//    runTime.writeAndEnd();

    // get the current time to reset time after mesh update
    scalar time = runTime.time().value();
    
    for (int i=1; i <= refineInterval; i++)
    {
//        for (int j=1; j <= maxRefinement; j++)
//        {
            if (!overwrite)
            {
                runTime++;
            }

            mesh.update();
//        }
    }

//    if(args.optionFound("overwrite"))
//    {
//        // mesh.update() calls write member if updateTime 
//        // if current runTime > first writeTime (of system/controlDict)
//        scalar latestTime = runTime.times().last().value();

//        // if write was called in mesh.update() delete the new time folder
//        if(latestTime > time)
//        {
//            system("rm -rf "+  runTime.timePath());
//        }

//        // reset time 
//        runTime.setTime(time, 0);
//    }

    runTime.writeNow();


    Info<< "End\n" << endl;

    mesh.clear();
    return(0);
}


// ************************************************************************* //
