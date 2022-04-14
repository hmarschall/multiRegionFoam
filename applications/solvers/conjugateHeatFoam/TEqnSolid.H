fvScalarMatrix* TsEqn = new fvScalarMatrix
(
    fvm::ddt(rhoS*cvS, Ts)
  ==
    fvm::laplacian(kS, Ts, "laplacian(k,Ts)")
);
