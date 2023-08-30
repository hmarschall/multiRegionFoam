// fuelCellFoam
Info<< "Constructing momentum equations" << endl;

fvVectorMatrix UEqn(U, rho.dimensions()*U.dimensions()*dimVol/dimTime);

// {
//     {
        UEqn =
        (
            phase_.UEqn()
        //  ==
        //     fvOptions(rho, U)
        );
        this->porosityModels().addResistance(UEqn);
        UEqn.relax();
        // fvOptions.constrain(UEqn);
        U.correctBoundaryConditions();
        // fvOptions.correct(U);
//     }
// }

Info << "Nach UEqn " << endl;

// // rhoPimpleFoam

// // Solve the Momentum equation

// tmp<fvVectorMatrix> UEqn
// (
//     phase_.UEqn()
// );

// this->porosityModels().addResistance(UEqn());

// UEqn().relax
// (
//     mesh.solutionDict().equationRelaxationFactor(U.select(pimple_.finalIter()))
// );

// volScalarField rUA = 1.0/UEqn().A();

// if (pimple_.momentumPredictor())
// {
//     Foam::solve
//     (
//         UEqn() == -fvc::grad(p_rgh_),
//         mesh.solutionDict().solver((U.select(pimple_.finalIter())))
//     );

//     // tmp<fvVectorMatrix> tUEqn
//     // (
//     //     UEqn() == -fvc::grad(p)
//     // );

//     // tUEqn.solve
//     // (
//     //     mesh.solutionDict().solver((U.select(pimple_.finalIter())))   
//     // );

// }
// else
// {
//     U = rUA*(UEqn().H() - fvc::grad(p_rgh_));
//     U.correctBoundaryConditions();
// }


// // ConstitutiveFoam


// fvVectorMatrix UEqn
// (
//     phase_.UEqn()
// );

// UEqn.relax();

// if (pimple_.momentumPredictor())
// {
//     Foam::solve(UEqn == -fvc::grad(p_rgh_));
// }
