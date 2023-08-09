{
    // Initialize MRF
    if (mesh().time().timeIndex() == 1)
    {
        if (movingReferenceFrame_.found("center"))
        {
            center_ = vector
            (
                movingReferenceFrame_.lookup("center")
            );

            center0_ = center_;
        }
        else
        {
            const fvMesh& mrfMesh = mesh().time().lookupObject<fvMesh>(centerFromMesh_);

            center_ = gSum(mrfMesh.C().internalField()*mrfMesh.V())/gSum(mrfMesh.V());

            movingReferenceFrame_.add("center", center_);

            center0_ = center_;
        }

        if (movingReferenceFrame_.found("XF"))
        {
            XF_ = dimensionedVector
            (
                movingReferenceFrame_.lookup("XF")
            );
        }
        else
        {
            XF_ = dimensionedVector
            (
                "XF",
                dimLength,
                center_
            );

            movingReferenceFrame_.add("XF", XF_);
        }

        if (movingReferenceFrame_.found("UF"))
        {
            UF_ = dimensionedVector
            (
                movingReferenceFrame_.lookup("UF")
            );
        }
        else
        {
            movingReferenceFrame_.add("UF", UF_);
        }

        if (movingReferenceFrame_.found("aF"))
        {
            aF_ = dimensionedVector
            (
                movingReferenceFrame_.lookup("aF")
            );
        }
        else
        {
            movingReferenceFrame_.add("aF", aF_);
        }

    }


    vector centerOld = center_;

    const fvMesh& mrfMesh = mesh().time().lookupObject<fvMesh>(centerFromMesh_);

    center_ = gSum(mrfMesh.C().internalField()*mrfMesh.V())/gSum(mrfMesh.V());

    Info << "\nCurrent MRF centre: " << center_ << endl;

    dimensionedVector dUF
    (
        "dUF",
        dimVelocity,
        lambdaFf_*(center_-center0_)/runTime().deltaT0().value()
      + lambdaF0_*(center0_-centerOld)/runTime().deltaT0().value()
    );

    XF_ += (UF_ + 0.5*dUF)*runTime().deltaT0();
    Info << "Current MRF position: " << XF_ << endl;

    UF_ += dUF;
    Info << "Current MRF velocity: " << UF_ << endl;

    aF_ = dUF/runTime().deltaT0();
    Info << "Current MRF acceleration: " << aF_ << endl;

    if (!closedVolume_ && hasSpacePatch_)
    {
        // get space patch index
        label spacePatchID = mesh().boundaryMesh().findPatchID("space");

        if
        (
            U_().boundaryField()[spacePatchID].type()
         == inletOutletFvPatchVectorField::typeName
        )
        {
            inletOutletFvPatchVectorField& spaceU =
                refCast<inletOutletFvPatchVectorField>
                (
                    U_().boundaryField()[spacePatchID]
                );

            spaceU.refValue() = -UF_.value();
        }
        else
        {
            FatalErrorIn("Foam::regionTypes::interTrackFluid::correct()")
                << "Velocity boundary condition at space patch is not "
                    << inletOutletFvPatchVectorField::typeName
                    << abort(FatalError);
        }
    }

    // Update movingReferenceFrame_ dictionary
    if (movingReferenceFrame_.found("center_"))
    {
        movingReferenceFrame_.remove("center_");
    }
    movingReferenceFrame_.add("center_", center_);

    if (movingReferenceFrame_.found("XF"))
    {
        movingReferenceFrame_.remove("XF");
    }
    movingReferenceFrame_.add("XF", XF_);

    if (movingReferenceFrame_.found("UF"))
    {
        movingReferenceFrame_.remove("UF");
    }
    movingReferenceFrame_.add("UF", UF_);

    if (movingReferenceFrame_.found("aF"))
    {
        movingReferenceFrame_.remove("aF");
    }
    movingReferenceFrame_.add("aF", aF_);
}