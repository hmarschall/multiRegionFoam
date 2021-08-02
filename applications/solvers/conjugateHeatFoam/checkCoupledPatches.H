// **********************************************************************
//  Lookup coupled patches
// **********************************************************************

    // Fluid
    forAll(T.boundaryField(), patchI)
    {
        if
        (
            T.boundaryField()[patchI].type() 
         == regionCoupleTemperatureFvPatchScalarField::typeName
        )
        {
            fluidCoupledTemperaturePatchIDs.append(patchI);
            Info << "Added patch " 
                 << T.mesh().boundaryMesh()[patchI].name() 
                 << " to fluidCoupledTemperaturePatchIDs."
                 << endl;
        }
        else if
        (
            T.boundaryField()[patchI].type() 
         == regionCoupleHeatFluxFvPatchScalarField::typeName
        )
        {
            fluidCoupledHeatFluxPatchIDs.append(patchI);
            Info << "Added patch " 
                 << T.mesh().boundaryMesh()[patchI].name() 
                 << " to fluidCoupledHeatFluxPatchIDs."
                 << endl;
        }
    }
    fluidCoupledTemperaturePatchIDs.shrink();
    fluidCoupledHeatFluxPatchIDs.shrink();

    // Solid
    forAll(Ts.boundaryField(), patchI)
    {
        if
        (
            Ts.boundaryField()[patchI].type() 
         == regionCoupleTemperatureFvPatchScalarField::typeName
        )
        {
            solidCoupledTemperaturePatchIDs.append(patchI);
            Info << "Added patch " 
                 << T.mesh().boundaryMesh()[patchI].name() 
                 << " to solidCoupledTemperaturePatchIDs."
                 << endl;
        }
        else if
        (
            Ts.boundaryField()[patchI].type() 
         == regionCoupleHeatFluxFvPatchScalarField::typeName
        )
        {
            solidCoupledHeatFluxPatchIDs.append(patchI);
            Info << "Added patch " 
                 << T.mesh().boundaryMesh()[patchI].name() 
                 << " to solidCoupledHeatFluxPatchIDs."
                 << nl << endl;
        }
    }
    solidCoupledTemperaturePatchIDs.shrink();
    solidCoupledHeatFluxPatchIDs.shrink();

//    Info << "Coupled patch IDs:" << nl
//         << " solid: "
//         << " temperature: " << solidCoupledTemperaturePatchIDs
//         << " heat flux: " << solidCoupledHeatFluxPatchIDs << nl
//         << " fluid: "
//         << " temperature: " << fluidCoupledTemperaturePatchIDs
//         << " heat flux: " << fluidCoupledHeatFluxPatchIDs << nl
//         << endl;


// **********************************************************************
//  Check coupled patches
// **********************************************************************

    int nCoupledTemperaturePatches =
        fluidCoupledTemperaturePatchIDs.size() +
        solidCoupledTemperaturePatchIDs.size();

    int nCoupledHeatFluxPatches =
        fluidCoupledHeatFluxPatchIDs.size() +
        solidCoupledHeatFluxPatchIDs.size();

    if((nCoupledTemperaturePatches != nCoupledHeatFluxPatches))
    {
        Info<< "\nError: Total number of regionCoupleTemperature patches "
            << "!= total number of regionCoupleHeatFlux patches "
            << nl << endl;
        exit(1);
    }

