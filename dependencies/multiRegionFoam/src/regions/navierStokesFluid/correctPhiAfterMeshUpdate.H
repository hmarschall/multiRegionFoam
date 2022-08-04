{
    for (label i = 0; i < p_.boundaryField().size(); i++)
    {
        if (p_.boundaryField()[i].fixesValue())
        {
            pcorrTypes_[i] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorrMesh",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("pcorrMesh", p_.dimensions(), 0.0),
        pcorrTypes_
    );

    dimensionedScalar rAUf
    (
        "(1|A(U))",
        dimTime/rho_.dimensions(),
        mesh().time().deltaT().value()
    );

    phi_ = (fvc::interpolate(U_) & mesh().Sf());

    //adjustPhi(phi_, U_, pcorr);

    mesh().schemesDict().setFluxRequired(pcorr.name());

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::div(phi_)
        );

        if (whichProcHasRef[Pstream::myProcNo()])
        {
            pcorrEqn.setReference(pRefCell_, pRefValue_);
        }

        pcorrEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi_ -= pcorrEqn.flux();
        }
    }
}
