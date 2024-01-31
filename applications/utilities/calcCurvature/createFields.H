// Create Finite Area mesh
faMesh aMesh(mesh);

// Detect the free surface or interface patch
label fsPatchID = -1;
forAll (mesh.boundary(), patchI)
{
    if(mesh.boundary()[patchI].name() == "freeSurface")
    {
	    fsPatchID = patchI;

		Info << "Found free surface patch. ID: " << fsPatchID
			 << endl;
	}
    else if(mesh.boundary()[patchI].name() == "interface")
    {
	    fsPatchID = patchI;

		Info << "Found interface patch. ID: " << fsPatchID
			 << endl;
	}
}

if(fsPatchID == -1)
{
    FatalErrorIn("createFields")
		<< abort(FatalError);
}

// Create the curvature field, curvatures only on the faMesh
volScalarField curvature
(
    IOobject
    (
        "curvature",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("curvature", dimless/dimLength, 0)
);

// Create the curvature field, curvatures only on the faMesh
volScalarField curvatureNew
(
    IOobject
    (
        "curvatureNew",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("curvatureNew", dimless/dimLength, 0)
);
