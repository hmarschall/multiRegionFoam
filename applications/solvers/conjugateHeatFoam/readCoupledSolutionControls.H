    IOdictionary coupledSolutionDict
    (
        IOobject
        (
            "coupledSolutionDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    couplingMethods couplingMethod = partitioned;

    word couplingMethodName = 
        coupledSolutionDict.lookupOrDefault<word>("couplingMethod", "partitioned");

    if (couplingMethodName == "partitioned")
    {
        couplingMethod = partitioned;
    }
    else if (couplingMethodName == "monolithic")
    {
        couplingMethod = monolithic;
    }
    else
    {
        FatalErrorIn("main(int argc, char *argv[])")
            << couplingMethodName << " not supported"
            << exit(FatalError);

        FatalError.exit(1);
    }

    scalar relConvergeTolerance =
        coupledSolutionDict.subDict("partitioned")
        .lookupOrDefault<scalar>("relConvergeTolerance", 0);

    int maxCoupleIter =
        coupledSolutionDict.subDict("partitioned")
        .lookupOrDefault<int>("maxCoupleIter", 50);
