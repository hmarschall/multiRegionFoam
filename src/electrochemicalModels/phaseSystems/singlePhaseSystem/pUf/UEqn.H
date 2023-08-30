Info<< "Constructing face momentum equations" << endl;

fvVectorMatrix UEqn(U, rho.dimensions()*U.dimensions()*dimVol/dimTime);

{
    {
        UEqn =
        (
            phase_.UfEqn()
        //  ==
        //     fvOptions(rho, U)
        );
        this->porosityModels().addResistance(UEqn);
        UEqn.relax();
        // fvOptions.constrain(UEqn);
        U.correctBoundaryConditions();
        // fvOptions.correct(U);
    }
}

