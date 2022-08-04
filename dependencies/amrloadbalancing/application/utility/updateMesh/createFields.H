     // read all registered field objects
    IOobjectList objects(mesh, runTime.timeName());

    // Read volScalarFields
    PtrList<volScalarField> vsFlds;
    ReadFields(mesh, objects, vsFlds);

    // Read volVectorFields
    PtrList<volVectorField> vvFlds;
    ReadFields(mesh, objects, vvFlds);

    // Read surfaceScalarFields (phi)
    PtrList<surfaceScalarField> ssFlds;
    ReadFields(mesh, objects, ssFlds);

    // Read surfaceVectorFields (Uf)
    PtrList<surfaceVectorField> svFlds;
    ReadFields(mesh, objects, svFlds);

    
    IOdictionary refineDict
    (
        IOobject
        (
            "dynamicMeshDict",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

//    word dynamicFvMesh; 
//    refineDict.lookup("dynamicFvMesh") >> dynamicFvMesh;

//    scalar refineInterval = readScalar
//    (
//        refineDict.subDict(dynamicFvMesh + "Coeffs")
//        .lookup("refineInterval")
//    );

    scalar refineInterval = readScalar
    (
        refineDict.subDict("dynamicRefineFvMeshCoeffs")
        .lookup("refineInterval")
    );

//    scalar maxRefinement = 
//        readScalar
//        (
//            refineDict.subDict("dynamicRefineFvMeshCoeffs")
//            .lookup("maxRefinement")
//        );
