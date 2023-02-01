{
    // Initialise flux with interpolated velocity
    phi_() = fvc::interpolate(U_()) & mesh().Sf();

    // //- CH: why is this needed?
    // if (mesh().changing())
    // {
    //     forAll(U_().boundaryField(), patchi)
    //     {
    //         if (U_().boundaryField()[patchi].fixesValue())
    //         {
    //             U_().boundaryField()[patchi].initEvaluate();
    //         }
    //     }

    //     forAll(U_().boundaryField(), patchi)
    //     {
    //         if (U_().boundaryField()[patchi].fixesValue())
    //         {
    //             U_().correctBoundaryConditions();
    //             U_().boundaryField()[patchi].evaluate();

    //             phi_().boundaryField()[patchi] =
    //                 U_().boundaryField()[patchi] 
    //                 & mesh().Sf().boundaryField()[patchi];
    //         }
    //     }
    // }

    wordList pcorrTypes(p_().boundaryField().types());

    for (label i=0; i<p_().boundaryField().size(); i++)
    {
        if (p_().boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("pcorr", p_().dimensions(), 0.0),
        pcorrTypes
    );

//     //CH: do we need to adjust phi here?
//     if (pcorr.needReference())
//     {
//         // Make flux relative
//         phi -= fvc::meshPhi(rho, U);
//         // Adjust/Scale phi
// #       include correctPatchPhi.H
//         // Make flux absolute
//         phi += fvc::meshPhi(rho, U);
//     }

    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian
            (
                rAU_,
                pcorr,
                "laplacian((1|A(U)),p)" 
            ) == fvc::div(phi_())
        );

        pcorrEqn.setReference(pRefCell_, pRefValue_);
        pcorrEqn.solve();

        if (pimple_.finalNonOrthogonalIter())
        {
            phi_() -= pcorrEqn.flux();
        }
    }
}

//#include "continuityErrs.H"
