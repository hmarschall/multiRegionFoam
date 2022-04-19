for (int Ecorr=0; Ecorr<nEnergyCorrectors; Ecorr++)
{
    phaseModel& phase = electrochemicalPhaseModel_();

    const volScalarField& alpha = phase;

    volVectorField& U = phase.URef();

    basicRhoThermo& thermo = phase.thermoRef();

    volScalarField& rho = thermo.rho();

    if (!phase.isothermal())
    {
        fvScalarMatrix EEqn
        (
            phase.heEqn()
        );

        EEqn.relax();
        // fvOptions.constrain(EEqn);
        EEqn.solve();
        // fvOptions.correct(thermo.h());
    }
}
