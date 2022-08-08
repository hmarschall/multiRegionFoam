fvScalarMatrix* TEqn = new fvScalarMatrix
(
    fvm::ddt(T)
  + fvm::div(phi, T)
  ==
    fvm::laplacian(alpha, T)
);
