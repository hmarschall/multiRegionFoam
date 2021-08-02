{
    // Initialize subiteration residuals
    scalar patchMaxTemperatureResidual = 0;
    scalar patchMaxHeatFluxResidual = 0;
    scalar globalMaxTemperatureResidual = 0;
    scalar globalMaxHeatFluxResidual = 0;


// **********************************************************************
//  Temperature coupling
// **********************************************************************

    // Fluid
    forAll(fluidCoupledTemperaturePatchIDs, patchI)
    {
        label patchID = fluidCoupledTemperaturePatchIDs[patchI];

        const fvPatchField<scalar>& coupledTemperaturePatchField = 
            T.boundaryField()[patchID];

        // Maximum normalized coupled patch temperature residual
        patchMaxTemperatureResidual = 
            refCast<const regionCoupleTemperatureFvPatchScalarField>
            (coupledTemperaturePatchField).maxResidual();

        globalMaxTemperatureResidual = 
            max
            (
                globalMaxTemperatureResidual,
                patchMaxTemperatureResidual
            );

        // Print patch total heat transfer rate
        scalar patchHeatRate =
            gSum
            (
                k.value() *
                coupledTemperaturePatchField.snGrad() *
                coupledTemperaturePatchField.patch().magSf()
            );

        Info<< "Total heat transfer rate (Fluid, patch = "
            << coupledTemperaturePatchField.patch().name() << "): "
            << patchHeatRate << " Watts" << endl;
    }

    // Solid
    forAll(solidCoupledTemperaturePatchIDs, patchI)
    {
        label patchID = solidCoupledTemperaturePatchIDs[patchI];

        const fvPatchField<scalar>& coupledTemperaturePatchField = 
            Ts.boundaryField()[patchID];

        // Maximum normalized coupled patch temperature residual
        patchMaxTemperatureResidual =
            refCast<const regionCoupleTemperatureFvPatchScalarField>
            (coupledTemperaturePatchField).maxResidual();

        globalMaxTemperatureResidual = 
            max
            (
                globalMaxTemperatureResidual,
                patchMaxTemperatureResidual
            );

        // Print patch total heat transfer rate
        scalar patchHeatRate =
            gSum
            (
                ksolid.value() *
                coupledTemperaturePatchField.snGrad() *
                coupledTemperaturePatchField.patch().magSf()
            );

        Info<< "Total heat transfer rate (Solid, patch = "
            << coupledTemperaturePatchField.patch().name() << "): "
            << patchHeatRate << " Watts" << endl;
    }


// **********************************************************************
//  Heat flux coupling
// **********************************************************************

    // Fluid
    forAll(fluidCoupledHeatFluxPatchIDs, patchI)
    {
        label patchID = fluidCoupledHeatFluxPatchIDs[patchI];

        const fvPatchField<scalar>& coupledHeatFluxPatchField = 
            T.boundaryField()[patchID];

        // Maximum normalized coupled patch heat flux residual
        patchMaxHeatFluxResidual =
            refCast<const regionCoupleHeatFluxFvPatchScalarField>
            (coupledHeatFluxPatchField).maxResidual();

        globalMaxHeatFluxResidual = 
            max
            (
                globalMaxHeatFluxResidual,
                patchMaxHeatFluxResidual
            );

        // Print patch total heat transfer rate
        scalar patchHeatRate =
            gSum
            (
                k.value() *
                coupledHeatFluxPatchField.snGrad() *
                coupledHeatFluxPatchField.patch().magSf()
            );

        Info<< "Total heat transfer rate (Fluid, patch = "
            << coupledHeatFluxPatchField.patch().name() << "): "
            << patchHeatRate << " Watts" << endl;
    }

    // Solid
    forAll(solidCoupledHeatFluxPatchIDs, patchI)
    {
        label patchID = solidCoupledHeatFluxPatchIDs[patchI];

        const fvPatchField<scalar>& coupledHeatFluxPatchField = 
            Ts.boundaryField()[patchID];

        // Maximum normalized coupled patch heat flux residual
        patchMaxHeatFluxResidual =
            refCast<const regionCoupleHeatFluxFvPatchScalarField>
            (coupledHeatFluxPatchField).maxResidual();

        globalMaxHeatFluxResidual = 
            max
            (
                globalMaxHeatFluxResidual,
                patchMaxHeatFluxResidual
            );

        // Print patch total heat transfer rate
        scalar patchHeatRate =
            gSum
            (
                ksolid.value() *
                coupledHeatFluxPatchField.snGrad() *
                coupledHeatFluxPatchField.patch().magSf()
            );

        Info<< "Total heat transfer rate (Solid, patch = "
            << coupledHeatFluxPatchField.patch().name() << "): "
            << patchHeatRate << " Watts" << endl;
    }

    Info<< "Coupled solution residuals:" << nl
        << "\tMaximum normalized temperature residual = "
        << globalMaxTemperatureResidual << nl
        << "\tMaximum normalized heat flux residual = "
        << globalMaxHeatFluxResidual << endl;

    if
    (
        (globalMaxTemperatureResidual < relConvergeTolerance) &&
        (globalMaxHeatFluxResidual < relConvergeTolerance)
    )
    {
        Info<< "Coupled solution convergence reached after "
            << coupleIter << " subiterations" << nl << endl;
        break;
    }
    else if(coupleIter == maxCoupleIter)
    {
        Info<< "Coupled solution failed to converge after "
            << coupleIter << " subiterations" << nl << endl;
        exit(1);
    }
}

