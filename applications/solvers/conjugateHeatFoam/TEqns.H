// Subiteration loop (used for multi-region coupling)
if (couplingMethod == partitioned)
{
    for (int coupleIter=1; coupleIter<=maxCoupleIter; coupleIter++)
    {
        Info<< nl << "Couple Iteration = " << coupleIter << endl;

        // Solve thermal energy transport in the "fluid"
        #include "readFluidControls.H"

        for (int nonOrth=0; nonOrth<=nNonOrthCorrFluid; nonOrth++)
        {
            #include "TEqnFluid.H"

            TEqn->relax();
            TEqn->solve();
        }

        // Solve thermal energy transport in the "solid"
        #include "readSolidControls.H"

        for (int nonOrth=0; nonOrth<=nNonOrthCorrSolid; nonOrth++)
        {
            #include "TEqnSolid.H"

            TsEqn->relax();
            TsEqn->solve();
        }

        // Check coupled convergence
        #include "convergenceCheck.H"
    }
}
else if (couplingMethod == monolithic)
{
    coupledFvScalarMatrix TEqns(2);

    #include "TEqnFluid.H"

    #include "TEqnSolid.H"

    // Add fluid equation
    TEqns.set(0, TEqn);

    // Add solid equation
    TEqns.set(1, TsEqn);

    TEqns.solve(mesh.solutionDict().solver(Ts.name() + "coupled"));
}
//else
//{
//     FatalError << "Coupling method either partitioned or monolithic."
//                << exit(FatalError);
//}

