
    fvVectorMatrix ddtUEqn(fvm::ddt(rho_(), U_()));

    mrfZones_.translationalMRFs().addFrameAcceleration(ddtUEqn, rho_());
    
    // Momentum equation
    fvVectorMatrix UEqn
    (
        ddtUEqn
      + fvm::div(fvc::interpolate(rho_())*phi_(), U_(), "div(phi,U)")
      - fvm::laplacian(mu_(), U_())
    );

    rAU_() = 1.0/UEqn.A();

    // Under-relax momentum. Note this will destroy the H and A
    UEqn.relax();

    // Insert momentum equation
    UpEqn.insertEquation(0, UEqn);

#   include "addBlockCoupledBC.H"

