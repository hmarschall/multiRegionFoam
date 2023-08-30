Info << "\nReading g" << endl;
const uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

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


volVectorField& U = phase_.URef();
surfaceScalarField& phi = phase_.phiRef();
// const surfaceScalarField& alphaPhi = phase_.alphaPhi();
surfaceScalarField& rhoPhi = phase_.alphaRhoPhiRef();

rhoThermo& thermo = phase_.thermoRef();

volScalarField& p = thermo.p();
volScalarField& rho = thermo.rho();
const volScalarField& psi = thermo.psi();
volScalarField dmdt = phase_.dmdt();

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    // *this
    readScalar(this->lookup("pMin"))
);

// fv::options& fvOptions = this->fvOptions();

// TODO: For fluxRequired Bug -> Important line!
mesh.schemesDict().setFluxRequired(p_rgh_.name());