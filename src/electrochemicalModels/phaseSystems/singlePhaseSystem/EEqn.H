int nEnergyCorrectors
(
    pimple_.dict().lookupOrDefault<int>("nEnergyCorrectors", 1)
);

for (int Ecorr=0; Ecorr<nEnergyCorrectors; Ecorr++)
{

    Info << "Solving in fluid region for " << phase_.name() << " EEqn " << endl;

    const volScalarField& alpha = phase_;

    if (!phase_.isothermal())
    {
        fvScalarMatrix EEqn
        (
            phase_.heEqn()
        //  ==
        //    alpha*rho*(U&g)
        );

        EEqn.relax();
        // fvOptions.constrain(EEqn);
        EEqn.solve();
        // fvOptions.correct(thermo.he());
    }
}
