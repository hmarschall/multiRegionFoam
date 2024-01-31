// Create Finite Area mesh
faMesh aMesh(mesh);

edgeScalarField nonOrtho
(
    IOobject
    (
        "nonOrtho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    aMesh,
    dimensionedScalar("zero", dimless, 0.0)
);

edgeScalarField skew
(
    IOobject
    (
        "skew",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    aMesh,
    dimensionedScalar("zero", dimless, 0.0)
);

edgeVectorField edgeSv
(
    IOobject
    (
        "edgeSv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    aMesh,
    dimensionedVector("zero", dimless, vector::zero)
);

edgeVectorField curvSv
(
    IOobject
    (
        "curvSv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    aMesh,
    dimensionedVector("zero", dimless, vector::zero)
);
