{
    volScalarField pcorr("pcorr", p_());
    pcorr *= 0;

    // Initialise flux with interpolated velocity
    phi_() = fvc::interpolate(U_()) & mesh().Sf();

    adjustPhi(phi_(), U_(), pcorr);

    mesh().schemesDict().setFluxRequired(pcorr.name());

    // while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAU_(), pcorr) == fvc::div(phi_())
        );

        pcorrEqn.setReference(pRefCell_, pRefValue_);
        pcorrEqn.solve();

        // if (pimple.finalNonOrthogonalIter())
        {
            phi_() -= pcorrEqn.flux();
        }

        // Fluxes are corrected to absolute velocity and further corrected
        // later.  HJ, 6/Feb/2009
    }

//#   include "continuityErrs.H"
}

