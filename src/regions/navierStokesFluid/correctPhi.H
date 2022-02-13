{
//    volScalarField pcorr("pcorr", p_);
//    pcorr *= 0;

    if (mesh().changing())
    {
        forAll(U_.boundaryField(), patchi)
        {
            if (U_.boundaryField()[patchi].fixesValue())
            {
                U_.boundaryField()[patchi].initEvaluate();
            }
        }

        forAll(U_.boundaryField(), patchi)
        {
            if (U_.boundaryField()[patchi].fixesValue())
            {
                U_.correctBoundaryConditions();
                U_.boundaryField()[patchi].evaluate();

                phi_.boundaryField()[patchi] =
                    U_.boundaryField()[patchi] 
                    & mesh().Sf().boundaryField()[patchi];
            }
        }
    }

    dimensionedScalar rAUf
    (
        "(1|A(U))",
        dimTime/dimMass*dimVolume,
        mesh().time().deltaT().value()
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian
            (
                rAUf,
                pcorr_,
                "laplacian((1|A(U)),p)" 
            ) == fvc::div(phi_)
        );

        pcorrEqn.setReference(pRefCell_, pRefValue_);
        pcorrEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi_ -= pcorrEqn.flux();
        }
    }
}

//#include "continuityErrs.H"
