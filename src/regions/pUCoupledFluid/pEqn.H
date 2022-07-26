// Pressure parts of the continuity equation
presSource_.reset
(
    new surfaceScalarField
    (
      "presSource",
      fvc::interpolate(rAU_())*
      (fvc::interpolate(fvc::grad(p_())) & mesh().Sf())
    )
);

fvScalarMatrix pEqn
(
    - fvm::laplacian(rAU_(), p_())
  ==
    - fvc::div(presSource_())
);

fvScalarMatrices.set
(
    p_().name() + mesh().name() + "Eqn",
    new fvScalarMatrix(pEqn)
);

pEqn.setReference(pRefCell_, pRefValue_);

UpEqn.insertEquation(3, pEqn);
