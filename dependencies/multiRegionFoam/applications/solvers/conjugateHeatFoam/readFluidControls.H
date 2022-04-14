    const dictionary& fluidControlsDict = mesh.solutionDict().subDict("SIMPLE");

    int nNonOrthCorrFluid = 0;

    if (fluidControlsDict.found("nNonOrthogonalCorrectors"))
    {
        nNonOrthCorrFluid = readInt(fluidControlsDict.lookup("nNonOrthogonalCorrectors"));
    }
