Info<< nl << "Reading thermophysicalProperties" << endl;

autoPtr<rhoReactionThermo> thermo = rhoReactionThermo::New(mesh);

basicMultiComponentMixture& composition = thermo->composition();
PtrList<volScalarField>& Y = composition.Y();

word inertSpecie(thermo->defaultSpecie());

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    thermo->rho()
);


Info<< nl << "Ending simpleREactionFoam " << endl;
// const volScalarField& psi = thermo->psi();
// volScalarField& h = thermo->h();
// const volScalarField& T = thermo->T();
