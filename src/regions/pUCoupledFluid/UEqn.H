
    fvVectorMatrix ddtUEqn(fvm::ddt(U_()));

    mrfZones_.translationalMRFs().addFrameAcceleration(ddtUEqn);
    
    // Momentum equation
    fvVectorMatrix UEqn
    (
        ddtUEqn
      + fvm::div(phi_(), U_())
      - fvm::laplacian(nu_(), U_())
    );

    rAU_() = 1.0/UEqn.A();

    // Under-relax momentum. Note this will destroy the H and A
    UEqn.relax();

    // Insert momentum equation
    UpEqn.insertEquation(0, UEqn);

#   include "addBlockCoupledBC.H"

