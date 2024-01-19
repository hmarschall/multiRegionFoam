// NON_ORTHOGONALITY

// For all internal edges check that the d dot S product is positive

// Severe nonorthogonality threshold
scalar nonOrthThreshold = 70; //deg
const scalar severeNonorthogonalityThreshold =
    ::cos(nonOrthThreshold/180.0*mathematicalConstant::pi);

scalar minDDotS = GREAT;
scalar sumDDotS = 0;
label severeNonOrth = 0;
label errorNonOrth = 0;

forAll (eNei, edgeI)
{
    vector d = aC[eNei[edgeI]] - aC[eOwn[edgeI]];
    const vector& s = edgeL[edgeI];

    scalar dDotS = (d & s)/(mag(d)*mag(s) + VSMALL);
    
    nonOrtho[edgeI] = Foam::acos(dDotS-SMALL)/mathematicalConstant::pi*180.0;

    if (dDotS < severeNonorthogonalityThreshold)
    {
        if (dDotS > SMALL)
        {
            severeNonOrth++;
        }
        else
        {
            errorNonOrth++;
        }
    }

    if (dDotS < minDDotS)
    {
        minDDotS = dDotS;
    }

    sumDDotS += dDotS;
}


// SKEWNESS

scalar skewThreshold = 4;
scalar maxSkew = 0;
label nWarnSkew = 0;

forAll(eNei, edgeI)
{
    vector Cpe = edgeCtrs[edgeI] - aC[eOwn[edgeI]];
    vector d = aC[eNei[edgeI]] - aC[eOwn[edgeI]];

    // Skewness vector
    vector sv =
        Cpe - ((edgeL[edgeI] & Cpe)/((edgeL[edgeI] & d) + SMALL))*d;
    
    skew[edgeI] = mag(sv);
    vector eN = aMesh.points()[faEdges[edgeI].end()] 
              - aMesh.points()[faEdges[edgeI].start()];
    edgeSv[edgeI] = (sv & eN)*eN;
    curvSv[edgeI] = sv - edgeSv[edgeI];
    
    // UNSURE:
    vector svHat = sv/(mag(sv) + VSMALL);
    // Normalisation distance calculated as the approximate distance
    // from the face centre to the edge of the face in the direction of
    // the skewness
    scalar fd = 0.2*mag(d) + VSMALL;
    fd = max(fd, mag(svHat & (aMesh.points()[faEdges[edgeI].start()] - edgeCtrs[edgeI])));
    fd = max(fd, mag(svHat & (aMesh.points()[faEdges[edgeI].end()] - edgeCtrs[edgeI])));

    // Normalised skewness
    scalar skewness = mag(sv)/fd;

    // Check if the skewness vector is greater than the PN vector.
    // This does not cause trouble but is a good indication of a poor mesh.
    if (skewness > skewThreshold)
    {
        nWarnSkew++;
    }

    if(skewness > maxSkew)
    {
        maxSkew = skewness;
    }
}

//------         
// VTK output file with edge fields: nonOrthogonality and skewness
// only for faMeshes without boundaries        
//------
if ( aMesh.boundary().size() == 0 )
{
    OFstream outputFile("meshQualityFields.vtk");
    outputFile.precision(12);
    outputFile  << "# vtk DataFile Version 2.0" << endl
                << "Test to visualize pointVectorField" << endl
                << "ASCII" << endl
                << "DATASET UNSTRUCTURED_GRID" << endl
                << "POINTS" << " " << faPoints.size() << " " << "float" <<endl;
    forAll(faPoints, pI)
    {
        scalar xp = faPoints[pI].x();
        scalar yp = faPoints[pI].y();
        scalar zp = faPoints[pI].z();
        outputFile << xp << " " << yp << " " << zp << endl;
    }
    outputFile  << "CELLS " << faEdges.size() << " " << faEdges.size()*3 << endl;
    forAll(faEdges, eI)
    {
        scalar labelStart = faEdges[eI].start();
        scalar labelEnd = faEdges[eI].end();
        outputFile << "2 " << labelStart << " " << labelEnd << endl;
    }
    outputFile  << "CELL_TYPES " << faEdges.size() << endl;
    forAll(faEdges, eI)
    {
        outputFile << "3" << endl;
    }

    outputFile  << "CELL_DATA " << faEdges.size() << endl
                << "SCALARS nonOrtho float 1" << endl 
                << "LOOKUP_TABLE default" << endl;
    forAll(faEdges, eI)
    {
        scalar nOrtho = nonOrtho[eI];
        outputFile << nOrtho << endl;
    }

    outputFile  << "SCALARS skewness float 1" << endl 
                << "LOOKUP_TABLE default" << endl;
    forAll(faEdges, eI)
    {
        scalar skn = skew[eI];
        outputFile << skn << endl;
    }
    outputFile  << "SCALARS edgeSkewness float 1" << endl 
                << "LOOKUP_TABLE default" << endl;
    forAll(faEdges, eI)
    {
        scalar ske = mag(edgeSv[eI]);
        outputFile << ske << endl;
    }
    outputFile  << "SCALARS curvSkewness float 1" << endl 
                << "LOOKUP_TABLE default" << endl;
    forAll(faEdges, eI)
    {
        scalar skc = mag(curvSv[eI]);
        outputFile << skc << endl;
    }
}
