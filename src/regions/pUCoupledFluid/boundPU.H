{
    // Bound the pressure
    dimensionedScalar p1 = min(p_());
    dimensionedScalar p2 = max(p_());

    if (p1 < pMin_ || p2 > pMax_)
    {
        Info<< "p: " << p1.value() << " " << p2.value()
            << ".  Bounding." << endl;

        p_().max(pMin_);
        p_().min(pMax_);
        p_().correctBoundaryConditions();
    }

    // Bound the velocity
    volScalarField magU = mag(U_());
    dimensionedScalar U1 = max(magU);

    if (U1 > UMax_)
    {
        Info<< "U: " << U1.value() << ".  Bounding." << endl;

        volScalarField Ulimiter = pos(magU - UMax_)*UMax_/(magU + smallU_)
            + neg(magU - UMax_);
        Ulimiter.max(scalar(0));
        Ulimiter.min(scalar(1));

        U_() *= Ulimiter;
        U_().correctBoundaryConditions();
    }
}
