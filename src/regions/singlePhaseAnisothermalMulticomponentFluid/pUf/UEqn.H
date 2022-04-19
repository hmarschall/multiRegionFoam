Info<< "Constructing face momentum equations" << endl;

fvVectorMatrix UEqn(U, rho.dimensions()*U.dimensions()*dimVol/dimTime);

{
    {
        UEqn =
        (
            electrochemicalPhaseModel_.UfEqn()
        //  ==
        //     fvOptions(rho, U)
        );
        electrochemicalPhaseModel_->porosityModels().addResistance(UEqn);
        UEqn.relax();
        // fvOptions.constrain(UEqn);
        U.correctBoundaryConditions();
        // fvOptions.correct(U);
    }
}

