/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

Application
    calcCurvature

Author
    Chiara Pesci <pesci@mma.tu-darmstadt.de>
    All rights reserved.

Description
    Compute the curvature on a finiteArea mesh

\*---------------------------------------------------------------------------*/

#include "dynamicFvMesh.H"
#include "fvCFD.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"
    //argList::validOptions.insert("diameter","scalar");
    timeSelector::addOptions(true, false);

    argList args(argc, argv);

#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);

//#   include "createDynamicFvMesh.H"
#   include "createNamedDynamicFvMesh.H"

#   include "createFields.H"

    Info << endl;

    Info << "Curvature computation" << endl;

    Info << endl;

    curvature.boundaryField()[fsPatchID] =
        aMesh.faceCurvatures().internalField();
    curvature.write();

// Evaluate interface normal and get curvature
        vectorField nf = mesh.boundary()[fsPatchID].nf();

        vectorField n =
            aMesh.S()*fac::edgeIntegrate
            (
                aMesh.Le()*aMesh.edgeLengthCorrection()
            )().internalField();

        forAll (n, faceI)
        {
            if (mag(n[faceI])>SMALL)
            {
                n[faceI] /= mag(n[faceI]);
            }
            else
            {
                n[faceI] = nf[faceI];
            }
        }

	vectorField surfTensionForce =
            fac::edgeIntegrate
            (
                aMesh.Le()*aMesh.edgeLengthCorrection()
            )().internalField();

        scalarField normalSurfaceTensionForceDensity =
        (
            n & surfTensionForce
        );

	curvatureNew.boundaryField()[fsPatchID]  =
            normalSurfaceTensionForceDensity;
	curvatureNew.write();

    Info << "Min/Max curvature: "
         << min(curvature.boundaryField()[fsPatchID]) << " / "
         << max(curvature.boundaryField()[fsPatchID]) << endl;
    Info << "Avg curvature: "
         << sum(curvature.boundaryField()[fsPatchID])
            /curvature.boundaryField()[fsPatchID].size() << endl;

    Info << "Min/Max curvatureNew: "
         << min(curvatureNew.boundaryField()[fsPatchID]) << " / "
         << max(curvatureNew.boundaryField()[fsPatchID]) << endl;
    Info << "Avg curvatureNew: "
         << sum(curvatureNew.boundaryField()[fsPatchID])
            /curvatureNew.boundaryField()[fsPatchID].size() << endl;

    Info << endl;

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
