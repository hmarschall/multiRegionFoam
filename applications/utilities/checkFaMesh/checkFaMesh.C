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
    checkFaMesh

Description
    Check finiteArea mesh in terms of face non-orthogonality and skewness.
     
Author
    Chiara Pesci <pesci@mma.tu-darmstadt.de>
    All rights reserved.

\*---------------------------------------------------------------------------*/

#include "dynamicFvMesh.H"
#include "timeSelector.H"
#include "faMesh.H"
#include "faCFD.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "setRootCase.H"

#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);
    
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;
        
#       include "createMesh.H"
        
#       include "createFields.H"

        Info<< "Checking faMesh topology" << nl
            << "    Faces:         " << aMesh.patch().faceCentres().size() << nl
            << "    Edges:         " << aMesh.patch().nEdges() << nl
            << "    Points:        " << aMesh.patch().nPoints() << nl
            << "    Area [m^2]:    " << sum(aMesh.S()).value() << nl << endl;


        Info<< "Checking faMesh geometry..." << nl;
        
        // Mesh addressing and data
        const edgeList& faEdges = aMesh.edges();
        const pointField& faPoints = aMesh.points();
        const faceList& faFaces = aMesh.faces();
        // Internal edges normal (m in Tukovic 2012)
        const edgeVectorField& edgeL = aMesh.Le();
        const areaVectorField& aC = aMesh.areaCentres();
        const labelList& eOwn = aMesh.edgeOwner();
        const labelList& eNei = aMesh.edgeNeighbour();
        const edgeVectorField& edgeCtrs = aMesh.edgeCentres();
        
        // Face area
        scalar minFaArea = min(aMesh.S()).value();
        scalar maxFaArea = max(aMesh.S()).value();
        
        // Edge length
        scalarField magLe(faEdges.size(), 0.0);
        forAll (faEdges, edgeI)
        {
            vector p0 = faPoints[faEdges[edgeI].start()];
            vector p1 = faPoints[faEdges[edgeI].end()];
            magLe[edgeI] = mag(p1 - p0);
        }
        scalar minLe = min(magLe);
        scalar maxLe = max(magLe);
        
        // Face aspect ratio
        scalarField faceAR(faFaces.size(), 0.0);
        forAll ( faFaces, faceI )
        {
            labelList edgesFI = aMesh.patch().faceEdges()[faceI];
            scalarField edgeFILen(edgesFI.size(), 0.0);
            forAll (edgesFI, edgeI)
            {
                vector p0 = faPoints[faEdges[edgesFI[edgeI]].start()];
                vector p1 = faPoints[faEdges[edgesFI[edgeI]].end()];
                edgeFILen[edgeI] = mag(p1 - p0);
            }
            faceAR[faceI] = max(edgeFILen)/min(edgeFILen);
        }

#       include "primitiveFaMeshCheck.H"

        Info<< "    faMesh bounding box " << boundBox(aMesh.points()) << nl
            << "    Max aspect ratio = " << max(faceAR) << nl
            << "    Minumum face area = " << minFaArea << nl
            << "    Maximum face area = " << maxFaArea << nl
            << "    Minumum edge length = " << minLe << nl
            << "    Maximum edge length = " << maxLe << nl;
        label neiSize = eNei.size();
        if (neiSize > 0)
        {
            Info<< "    Mesh non-orthogonality Max: "
                << ::acos(minDDotS)/mathematicalConstant::pi*180.0
                << " average: " <<
                ::acos(sumDDotS/neiSize)/mathematicalConstant::pi*180.0
                << " Threshold = " << nonOrthThreshold
                << endl;
        }

        if (severeNonOrth > 0)
        {
            Info<< "   *Number of severely non-orthogonal faces: "
                << severeNonOrth << "." << endl;
        }            
        if (nWarnSkew > 0)
        {
            Info<< " ***Max skewness = " << maxSkew
                << ", " << nWarnSkew << " highly skew faces detected"
                << " Threshold = " << skewThreshold
                << endl;
        }
        else
        {
            Info<< "    Max skewness = " << maxSkew << " OK." << endl;
        }            
        Info<< "    Mesh skewness vector Max: " << max(skew.internalField()) << nl
            << "    Mesh skewness edge dir Max: " << max(mag(edgeSv.internalField())) << nl
            << "    Mesh skewness curvature Max: " << max(mag(curvSv.internalField())) << nl;
            //<< "    Planarity: " << tp << nl;
        
        Info<< endl;
        
        //------         
        // CSV output file with error statistics        
        //------ 
        // Open the error file for measurement output..
        OFstream errorFile("faMeshQuality.csv");
        // Nf : number of faces of the surface mesh
        // Ne : number of edges of the surface mesh
        // Np : number of points of the surface mesh
        // maxAR : maximum face aspect ratio
        // minAf : minimum face area
        // maxAf : maximum face area
        // minLe : minimum edge length
        // maxLe : maximum edge length
        // maxNonOrtho : maximum non-orthogonality in deg
        // avgNonOrtho : average non-orthogonality in deg
        // maxSkew : max skewness computed as checkMesh
        // maxSkVec : max skewness vector magnitude

        errorFile << "Nf,Ne,Np,maxAR,minAf,maxAf,minLe,maxLe,maxNonOrtho,avgNonOrtho,maxSkew,maxSkVec\n";

        errorFile << aMesh.patch().faceCentres().size() << ","
                << aMesh.patch().nEdges() << ","
                << aMesh.patch().nPoints() << ","
                << max(faceAR) << "," << minFaArea << "," << maxFaArea << "," 
                << minLe << "," << maxLe << ","
                << ::acos(minDDotS)/mathematicalConstant::pi*180.0 << ","
                << ::acos(sumDDotS/neiSize)/mathematicalConstant::pi*180.0 << ","
                << maxSkew << "," << max(skew.internalField()) << nl;
    
    } //-End time loop

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //

