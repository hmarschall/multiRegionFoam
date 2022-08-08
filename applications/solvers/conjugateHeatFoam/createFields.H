// **********************************************************************
//  Fluid side
// **********************************************************************

    //- Create fields
	Info<< "\tReading fluid temperature field\n" << endl;
	volScalarField T
	(
		IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	Info<< "\tReading fluid velocity field\n" << endl;
	volVectorField U
	(
		IOobject
		(
			"U",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh
	);

	Info<< "\tReading/calculating fluid face flux field\n" << endl;
	surfaceScalarField phi
	(
		IOobject
		(
			"phi",
			runTime.timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		linearInterpolate(U) & mesh.Sf()
	);


    //- transport properties
    IOdictionary fluidTransportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "\tReading fluid density\n" << endl;
    dimensionedScalar rho
    (
        fluidTransportProperties.lookup("rho")
    );

    Info<< "\tReading fluid specific heat\n" << endl;
    dimensionedScalar cp
    (
        fluidTransportProperties.lookup("cp")
    );

	Info<< "\tReading fluid thermal conductivity\n" << endl;
	dimensionedScalar k
	(
		fluidTransportProperties.lookup("k")
	);

    volScalarField* alphaPtr = NULL;

    if (couplingMethod == monolithic)
    {
        Info<< "Reading field alpha\n" << endl;

        alphaPtr = new volScalarField
        (
            IOobject
            (
                "alpha",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
    }
    else if (couplingMethod == partitioned)
    {
        alphaPtr = new volScalarField
        (
            IOobject
            (
                "alpha",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            k/(rho*cp)
        );
    }

    volScalarField& alpha = *alphaPtr;

    alpha = k/(rho*cp);

// **********************************************************************
//  Solid side
// **********************************************************************

    //- Create fields
	Info<< "\tReading solid temperature field\n" << endl;
    volScalarField Ts
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            solidMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        solidMesh
    );

    //- transport properties
    IOdictionary solidTransportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            solidMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

	Info<< "\tReading solid density\n" << endl;
	dimensionedScalar rhoS
	(
		solidTransportProperties.lookup("rho")
	);

	Info<< "\tReading solid specific heat\n" << endl;
	dimensionedScalar cvS
	(
		solidTransportProperties.lookup("cv")
	);

    Info<< "\tReading solid thermal conductivity\n" << endl;
    dimensionedScalar ksolid
    (
        solidTransportProperties.lookup("k")
    );

    volScalarField* kPtr = NULL;

    if (couplingMethod == monolithic)
    {
        Info<< "Reading field alpha\n" << endl;

        kPtr = new volScalarField
        (
            IOobject
            (
                "kS",
                runTime.timeName(),
                solidMesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidMesh
        );
    }
    else if (couplingMethod == partitioned)
    {
        kPtr = new volScalarField
        (
            IOobject
            (
                "alpha",
                runTime.timeName(),
                solidMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidMesh,
            ksolid
        );
    }

    volScalarField& kS = *kPtr;

    kS = ksolid;

