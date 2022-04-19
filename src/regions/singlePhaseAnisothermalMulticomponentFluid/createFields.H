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
tmp<volScalarField> mu = electrochemicalPhaseModel_->mu();
surfaceScalarField& phi = electrochemicalPhaseModel_->phiRef();
// const surfaceScalarField& alphaPhi = electrochemicalPhaseModel_.alphaPhi();
surfaceScalarField& rhoPhi = electrochemicalPhaseModel_->alphaRhoPhiRef();

basicRhoThermo& thermo = electrochemicalPhaseModel_->thermoRef();

volScalarField& p = thermo.p();
volScalarField& rho = thermo.rho();
const volScalarField& psi = thermo.psi();
volScalarField dmdt = electrochemicalPhaseModel_->dmdt();

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    // *this
    readScalar(dict_.lookup("pMin"))
);


// inverse matrix diagonal field rUA
volScalarField rUA
(
    IOobject
    (
        "(1|A(U))",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("rUA", dimTime/dimDensity, 1),
    zeroGradientFvPatchScalarField::typeName
);

// Pressure reference point in air
label pRefCell = 0;
scalar pRefValue = 0.0;

setRefCell
(
    p,
    mesh.solutionDict().subDict("PISO"),
    pRefCell,
    pRefValue
);


