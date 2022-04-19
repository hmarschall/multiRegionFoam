{
    if (!electrochemicalPhaseModel_->pure())
    {
        UPtrList<volScalarField>& Y = electrochemicalPhaseModel_->YActiveRef();

        forAll(Y, i)
        {
            fvScalarMatrix YiEqn
            (
                electrochemicalPhaseModel_->YiEqn(Y[i])
            //  ==
            //     fvOptions(rho, Y1[i])
            );

            YiEqn.relax();
            YiEqn.solve();

            Y[i].max(SMALL);
            Y[i].min(1.0);

            Info << "Y: " << Y[i].name() << ":"
                 << gMin(Y[i]) << " -> " << gMax(Y[i])
                 << endl;
        }
    }
}
