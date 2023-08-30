

    Info << "Solving in fluid region for " << phase_.name() << " TEqn " << endl;

    if (!phase_.isothermal())
    {
        fvScalarMatrix TEqn
        (
            phase_.TEqn()
        //  ==
        //    alpha*rho*(U&g)
        );

        TEqn.relax();
        // fvOptions.constrain(TEqn);
        TEqn.solve();
        // fvOptions.correct(thermo.he());
    }
