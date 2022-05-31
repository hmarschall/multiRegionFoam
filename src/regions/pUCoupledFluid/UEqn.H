
    // Momentum equation
    fvVectorMatrix UEqn
    (
        fvm::ddt(U_())
      + fvm::div(phi_(), U_())
      - fvm::laplacian(nu_(), U_())
    );

    rAU_() = 1.0/UEqn.A();

    // Under-relax momentum. Note this will destroy the H and A
    UEqn.relax();

    // Insert momentum equation
    UpEqn.insertEquation(0, UEqn);

#   include "addBlockCoupledBC.H"

