Info<< "Constructing face momentum equations" << endl;

fvVectorMatrix UEqn(U, rho.dimensions()*U.dimensions()*dimVol/dimTime);

{
    {
        UEqn =
        (
            electrochemicalPhaseModel_->UfEqn()
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


// Info << "Solving fluid momentum equations" << endl;

// p.storePrevIter();

// // Read controls
// const dictionary& piso = this->solutionDict().subDict("PISO");
// int nCorr(readInt(piso.lookup("nCorrectors")));
// int nNonOrthCorr(readInt(piso.lookup("nNonOrthogonalCorrectors")));

// fvVectorMatrix UEqn
// (
//     fvm::div(phi, U)
//     - fvm::laplacian(mu(), U)
// );

// electrochemicalPhaseModel_->porosityModels().addResistance(UEqn);

// UEqn.relax();

// solve(UEqn == -fvc::grad(p));
