        if (closedVolume_)
        {
            label intPatchID_ = 
                mesh().boundaryMesh().findPatchID("interfaceShadow");

            forAll(phi_().boundaryField(), patchI)
            {
                if
                (
                    !phi_().boundaryField()[patchI].coupled()
                 && patchI != intPatchID_
//                 && isA<zeroGradientFvPatchVectorField>
//                    (
//                        U_().boundaryField()[patchI]
//                    )
                )
                {
                    phi_().boundaryField()[patchI] ==
                    (
                        U_().boundaryField()[patchI]
                        & mesh().Sf().boundaryField()[patchI]
                    );
                }
            }
        }

        if (!closedVolume_)
        {
            forAll(phi_().boundaryField(), patchI)
            {
                if
                (
                   !phi_().boundaryField()[patchI].coupled()
                 && !p_().boundaryField()[patchI].fixesValue()
                )
                {
                    phi_().boundaryField()[patchI] ==
                    (
                        U_().boundaryField()[patchI]
                        & mesh().Sf().boundaryField()[patchI]
                    );
                }
            }
        }

        if (closedVolume_)
        {
            label intPatchID_ = 
                mesh().boundaryMesh().findPatchID("interfaceShadow");

            phi_().boundaryField()[intPatchID_] =
            (
                U_().boundaryField()[intPatchID_]
                & mesh().Sf().boundaryField()[intPatchID_]
            );

            scalarField weights =
                mag(phi_().boundaryField()[intPatchID_]);

            if(mag(gSum(weights)) > VSMALL)
            {
                weights /= gSum(weights);
            }

            phi_().boundaryField()[intPatchID_] -=
                weights*gSum(phi_().boundaryField()[intPatchID_]);

            phi_().boundaryField()[intPatchID_] +=
                p_().boundaryField()[intPatchID_].snGrad()
               *mesh().magSf().boundaryField()[intPatchID_]
               /AU.boundaryField()[intPatchID_];
        }

        if (!closedVolume_ && hasSpacePatch_)
        {
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
