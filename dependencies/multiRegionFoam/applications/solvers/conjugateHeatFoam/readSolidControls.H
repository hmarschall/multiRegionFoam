    const dictionary& solidControlsDict = solidMesh.solutionDict().subDict("SIMPLE");

    int nNonOrthCorrSolid = 0;
    if (solidControlsDict.found("nNonOrthogonalCorrectors"))
    {
        nNonOrthCorrSolid = readInt(solidControlsDict.lookup("nNonOrthogonalCorrectors"));
    }
