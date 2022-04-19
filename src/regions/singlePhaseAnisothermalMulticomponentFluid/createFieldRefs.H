const uniformDimensionedVectorField& g =
  runTime.lookupObject<uniformDimensionedVectorField>("g");


Info<< "\nReading hRef" << endl;
uniformDimensionedScalarField hRef
(
    IOobject
    (
        "hRef",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedScalar("hRef", dimLength, 0)
);

Info<< "Calculating field g.h\n" << endl;
dimensionedScalar ghRef
(
    mag(g.value()) > SMALL
    ? g & (cmptMag(g.value())/mag(g.value()))*hRef
    : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
);

volScalarField gh("gh", (g & mesh.C()) - ghRef);
surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);


volVectorField& U = electrochemicalPhaseModel_->URef();
surfaceScalarField& phi = electrochemicalPhaseModel_->phiRef();
// const surfaceScalarField& alphaPhi = electrochemicalPhaseModel_.alphaPhi();
surfaceScalarField& rhoPhi = electrochemicalPhaseModel_->alphaRhoPhiRef();

basicRhoThermo& thermo = electrochemicalPhaseModel_->thermoRef();

volScalarField& p = thermo.p();
volScalarField& rho = thermo.rho();
const volScalarField& psi = thermo.psi();
// volScalarField dmdt = electrochemicalPhaseModel_.dmdt();

// dimensionedScalar pMin
// (
//     "pMin",
//     dimPressure,
//     *this
// );
