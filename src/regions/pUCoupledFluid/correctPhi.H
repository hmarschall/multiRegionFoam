{
    pcorr_() *= 0;

    // Initialise flux with interpolated velocity
    phi_() = fvc::interpolate(U_()) & mesh().Sf();


    if (!closedVolume_ && hasSpacePatch_)
    {
        adjustPhi(phi_(), U_(), pcorr_());

        // get space patch index
        label scalePatchID =
            mesh().boundaryMesh().findPatchID("space");

        //- Non-const access to flux on patch
        fvsPatchField<scalar>& phip = 
            const_cast<fvsPatchField<scalar>& >
            (
                mesh().lookupObject<surfaceScalarField>("phi")
                .boundaryField()[scalePatchID]
            );

        scalar inletFlux = gSum(neg(phip)*phip);

        scalar outletFlux = gSum(pos(phip)*phip);

        if(outletFlux < VSMALL)
        {
            outletFlux = VSMALL;
        }

        scalar outflowScaling = -inletFlux/outletFlux;

        phip += pos(phip)*phip*(outflowScaling - 1.0);
    }

    mesh().schemesDict().setFluxRequired(pcorr_().name());

    // while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAU_(), pcorr_()) == fvc::div(phi_())
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

