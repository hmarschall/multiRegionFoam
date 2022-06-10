/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "freeSurfaceVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "fvCFD.H"
#include "faCFD.H"
#include "faMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::freeSurfaceVelocityFvPatchVectorField::
freeSurfaceVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    nonOrthCorr_(false),
    secondOrder_(false),
    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            this->db().time().constant(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    sigma_
    (
        transportProperties_.lookup("sigma")
    ),
    sigmaPtr_(NULL),
    curTimeIndex_(-1),
    myTimeIndex_( dimensionedInternalField().mesh().time().timeIndex() ),
    Fc_( p.patch().size(), vector::zero ),
    oldFc_( p.patch().size(), vector::zero ),
    oldoldFc_( p.patch().size(), vector::zero )
{}


Foam::freeSurfaceVelocityFvPatchVectorField::
freeSurfaceVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    nonOrthCorr_(false),
    secondOrder_(false),
    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            this->db().time().constant(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    sigma_
    (
        transportProperties_.lookup("sigma")
    ),
    sigmaPtr_(NULL),
    curTimeIndex_(-1),
    myTimeIndex_( dimensionedInternalField().mesh().time().timeIndex() ),
    Fc_( p.patch().size(), vector::zero ),
    oldFc_( p.patch().size(), vector::zero ),
    oldoldFc_( p.patch().size(), vector::zero )
{
    fvPatchVectorField::operator=(patchInternalField());

    Fc_ = p.patch().faceCentres();
    oldFc_ = p.patch().faceCentres();
    oldoldFc_ = p.patch().faceCentres();

    if (dict.found("secondOrder"))
    {
        secondOrder_ = Switch(dict.lookup("secondOrder"));
        Info << "Second order correction: " << secondOrder_ << endl;
    }

    if (dict.found("nonOrthCorr"))
    {
        nonOrthCorr_ = Switch(dict.lookup("nonOrthCorr"));
        Info << "NonOrthogonal correction: " << nonOrthCorr_ << endl;
    }
}


Foam::freeSurfaceVelocityFvPatchVectorField::
freeSurfaceVelocityFvPatchVectorField
(
    const freeSurfaceVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    nonOrthCorr_(ptf.nonOrthCorr_),
    secondOrder_(ptf.secondOrder_),
    transportProperties_(ptf.transportProperties_),
    sigma_(ptf.sigma_),
    sigmaPtr_(ptf.sigmaPtr_),
    curTimeIndex_(-1),
    myTimeIndex_( ptf.myTimeIndex_ ),
    Fc_( p.patch().size(), vector::zero ),
    oldFc_( p.patch().size(), vector::zero ),
    oldoldFc_( p.patch().size(), vector::zero )
{}


Foam::freeSurfaceVelocityFvPatchVectorField::
freeSurfaceVelocityFvPatchVectorField
(
    const freeSurfaceVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(fcvpvf, iF),
    phiName_(fcvpvf.phiName_),
    rhoName_(fcvpvf.rhoName_),
    nonOrthCorr_(fcvpvf.nonOrthCorr_),
    secondOrder_(fcvpvf.secondOrder_),
    transportProperties_(fcvpvf.transportProperties_),
    sigma_(fcvpvf.sigma_),
    sigmaPtr_(fcvpvf.sigmaPtr_),
    curTimeIndex_(-1),
    myTimeIndex_( fcvpvf.myTimeIndex_ ),
    Fc_( fcvpvf.oldFc_ ),
    oldFc_( fcvpvf.oldFc_ ),
    oldoldFc_( fcvpvf.oldoldFc_ )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::freeSurfaceVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

//    const fvMesh& mesh =
//        this->patch().boundaryMesh().mesh();

//    const volVectorField& U =
//        mesh.objectRegistry::lookupObject<volVectorField>("U");

//    label patchID = this->patch().index();


//    // Create Finite Area mesh
//    faMesh aMesh(mesh);

//    // Get surface areas
//    const scalarField& S = aMesh.S();

//    // Evaluate interface normal and get curvature
//    vectorField nf = this->patch().nf();

//    vectorField n = 
//        -S*fac::edgeIntegrate
//        (
//            aMesh.Le()*aMesh.edgeLengthCorrection()
//        )().internalField();

//    forAll (n, faceI)
//    {
//        if (mag(n[faceI])>SMALL)
//        {
//            n /= mag(n);
//        }
//        else
//        {
//            n[faceI] = nf[faceI];
//        }
//    }

//    // Correct normal velocity gradient at the free-surface
//    tmp<volVectorField> tUs
//    (
//        new volVectorField
//        (
//            IOobject
//            (
//                "Us",
//                db().time().timeName(),
//                mesh,
//                IOobject::NO_READ,
//                IOobject::NO_WRITE
//            ),
//            mesh,
//            vector::zero
//        )
//    );
//    volVectorField& Us = tUs();


//    string schemeSpec("skewCorrected");
//    tmp<fv::snGradScheme<vector> > tschemeV(
//        fv::snGradScheme<vector>::New(
//            U.mesh(),
//            IStringStream(schemeSpec)()
//        )
//    );
//    fv::snGradScheme<vector>& schemeV = tschemeV();

//    scalarField nSnGradU = (n & schemeV.snGrad(U)().boundaryField()[patchID]);

//    // Evaluate viscous and surface tension forces
//    IOdictionary transportProperties
//    (
//        IOobject
//        (
//            "surfaceProperties",
//            mesh.time().constant(),
//            mesh,
//            IOobject::MUST_READ,
//            IOobject::NO_WRITE
//        )
//    );
//    dimensionedScalar muFluidA(transportProperties.lookup("muFluidA"));

//    if 
//    (
//        mesh.foundObject<areaScalarField>("surfaceTension")
//    )
//    {
//        // contaminated interface
//        sigmaPtr_ = const_cast<areaScalarField*> 
//            (&mesh.lookupObject<areaScalarField>("surfaceTension"));
//    }
//    else
//    {
//        // clean interface
//        sigmaPtr_ = new areaScalarField
//        (
//            IOobject
//            (
//                "sigma",
//                this->db().time().timeName(),
//                aMesh.thisDb(),
//                IOobject::NO_READ,
//                IOobject::NO_WRITE
//            ),
//            aMesh,
//            sigma_,
//            zeroGradientFaPatchVectorField::typeName
//        );
//    }

//    areaScalarField& surfaceTension = *sigmaPtr_;

//    vectorField surfTensionForce(n.size(), vector::zero);

//    surfTensionForce =
//        S*fac::edgeIntegrate
//        (
//            linearEdgeInterpolate(surfaceTension)
//            *aMesh.Le()*aMesh.edgeLengthCorrection()
//        )().internalField();

//    vectorField tangentialSurfaceTensionForceDensity = 
//    (
//        (I-sqr(n)) & surfTensionForce
//    )/S;


//    // Set patch-normal gradient
//    gradient() = tangentialSurfaceTensionForceDensity/muFluidA.value()
//        + n*nSnGradU
//        - schemeV.snGrad(Us)().boundaryField()[patchID];

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    faMesh aMesh(mesh);

    const areaVectorField& nf = aMesh.faceAreaNormals();

    // Get interfacial curvature
    const areaScalarField& K = aMesh.faceCurvatures();

    // surface tension
    if
    (
        mesh.foundObject<areaScalarField>("surfaceTension")
    )
    {
        // contaminated interface
        sigmaPtr_ = const_cast<areaScalarField*>
            (&mesh.lookupObject<areaScalarField>("surfaceTension"));
    }
    else
    {
        // clean interface
        sigmaPtr_ = new areaScalarField
        (
            IOobject
            (
                "sigma",
                this->db().time().timeName(),
                aMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh,
            sigma_,
            zeroGradientFaPatchVectorField::typeName
        );
    }

    areaScalarField& sigma = *sigmaPtr_;

    areaVectorField gradSsigma =
    (
        fac::grad(sigma)// + aMesh.faceCurvatures()*sigma*nf - K*sigma*nf
    );

    // surface velocity terms
    areaVectorField Us
    (
        IOobject
        (
            "Us",
            db().time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh,
        vector::zero
    );

    const volVectorField& U =
        mesh.objectRegistry::lookupObject<volVectorField>("U");

    Us.internalField() = U.boundaryField()[this->patch().index()];

    areaScalarField divSU =
    (
        fac::div(Us)// + aMesh.faceCurvatures()*(nf&Us) - K*(nf&Us)
    );
    divSU.correctBoundaryConditions();

    areaTensorField gradSU =
    (
        fac::grad(Us)// + aMesh.faceCurvatures()*(Us*nf) - K*(Us*nf)
    );

    // Remove component of gradient normal to surface (area)
    gradSU -= nf*(nf & gradSU);
    gradSU.correctBoundaryConditions();

    dimensionedScalar mu(transportProperties_.lookup("mu"));

    // Set patch-normal gradient
    gradient() =
    (
        gradSsigma.internalField()/mu.value()
      + (
            divSU*nf
          + (gradSU&nf)
        )().internalField()
    );

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void Foam::freeSurfaceVelocityFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!updated())
    {
        updateCoeffs();
    }

//     fixedGradientFvPatchVectorField::evaluate();
    
    // Evaluate tangential component of velocity
    // using (second order discretisation) and nonorthogonal correction
    {
        vectorField n = patch().nf();
        vectorField delta = patch().delta();
        vectorField k = delta - n*(n&delta);

        word UName = this->dimensionedInternalField().name();

        tensorField gradUp(this->patch().size(), tensor::zero);

        if 
        (
            this->patch().boundaryMesh().mesh()
            .foundObject<volTensorField>("grad(" + UName + ")")
        )
        {
            const fvPatchField<tensor>& gradU =
                patch().lookupPatchField<volTensorField, tensor>
                (
                    "grad(" + UName + ")"
                );

            gradUp = gradU.patchInternalField();
        }
        else
        {
            const volVectorField& U = this->db().lookupObject<volVectorField>(UName);

            gradUp = patch().patchField<volTensorField, tensor>(fvc::grad(U));
//            gradUp = fvc::grad(U)().boundaryField()[this->patch().index()];
        }

        vectorField dUP(this->patch().size(), vector::zero);
        if (nonOrthCorr_)
        {
            dUP = (k&gradUp);
        }

        if (secondOrder_)
        {
            vectorField nGradUP = (n&gradUp);

            Field<vector>::operator=
            (
                this->patchInternalField() + dUP
              + 0.5*(gradient() + nGradUP)
               /this->patch().deltaCoeffs()
            );
        } 
        else
        {
            Field<vector>::operator=
            (
                this->patchInternalField() + dUP
              + gradient()/this->patch().deltaCoeffs()
            );
        }
    }


    // Calculate normal component of velocity from flux
    if (!this->db().objectRegistry::found(phiName_))
    {
        // Flux not available, do not update
        return;
    }

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

//    const fvsPatchField<scalar>& phip =
//        patch().patchField<surfaceScalarField, scalar>(phi);

//    scalarField phip = patch().patchField<surfaceScalarField, scalar>(phi);

//    vectorField n = patch().nf();
//    const Field<scalar>& magS = patch().magSf();

    const fvMesh & mesh = dimensionedInternalField().mesh();
    const fvPatch & p = patch();
    const polyPatch & pp = p.patch();

//    word ddtScheme
//    (
//        mesh.schemesDict().ddtScheme
//        (
//            "ddt(" + rho.name() + ',' + U().name()+')'
//        )
//    );

//    if
//    (
//        ddtScheme != fvbdf2DdtScheme<vector>::typeName
////     != fv::backwardDdtScheme<vector>::typeName
//    )

    if (myTimeIndex_ < mesh.time().timeIndex())
    {
        oldoldFc_ = oldFc_;
        oldFc_ = Fc_;
        myTimeIndex_ = mesh.time().timeIndex();
    }

    Fc_ = pp.faceCentres();

    // const pointField& oldPoints = mesh.oldPoints();
    const volVectorField & U = 
        mesh.lookupObject<volVectorField>(dimensionedInternalField().name());

    scalar deltaT = mesh.time().deltaT().value();
    scalar deltaT0 = mesh.time().deltaT0().value();

    if
    (
        U.oldTime().timeIndex() == U.oldTime().oldTime().timeIndex() 
     || U.oldTime().oldTime().timeIndex() < 0
    )
    {
        deltaT0 = GREAT;
    }

    // Set coefficients based on deltaT and deltaT0
    scalar coefft = 1 + deltaT / (deltaT + deltaT0);
    scalar coefft00 = deltaT * deltaT / ( deltaT0 * (deltaT + deltaT0) );
    scalar coefft0 = coefft + coefft00;

    tmp<vectorField> Up = (coefft*Fc_ - coefft0*oldFc_ + coefft00*oldoldFc_)
        /mesh.time().deltaT().value();

//    phip += p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U));
    scalarField phip = p.patchField<surfaceScalarField, scalar>
    (
        fvc::meshPhi(U)
    );

    tmp<vectorField> n = p.nf();
    const scalarField & magSf = p.magSf();
//    tmp<scalarField> Un = phip/(magSf + VSMALL);



    if (phi.dimensions() == dimVelocity*dimArea)
    {
        operator==(Up() + n()*(phip/(magSf + VSMALL) - (n()&Up())));
//        operator==(*this - n*(n & *this) + n*phip/magS);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        if (!this->db().objectRegistry::found(rhoName_))
        {
            // Rho not available, do not update
            return;
        }

        const fvPatchField<scalar>& rhop =
            lookupPatchField<volScalarField, scalar>(rhoName_);

        operator==(Up() + n()*(phip/(rhop*(magSf + VSMALL)) - (n()&Up())));
//        operator==(*this - n*(n & *this) + n*phip/(rhop*magS));
    }
    else
    {
        FatalErrorIn
        (
            "freeSurfaceVelocityFvPatchVectorField::evaluate()"
        )
            << "dimensions of phi are incorrect\n"
            << "    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

//    if (!this->db().objectRegistry::found(phiName_))
//    {
//        // Flux not available, do not update
//        return;
//    }

//    const surfaceScalarField& phi =
//        db().lookupObject<surfaceScalarField>(phiName_);

//    const fvsPatchField<scalar>& phip =
//        patch().patchField<surfaceScalarField, scalar>(phi);

//    const vectorField n(patch().nf());
//    const Field<scalar>& magS = patch().magSf();

//    if (phi.dimensions() == dimVelocity*dimArea)
//    {
//        operator==(*this - n*(n & *this) + n*phip/magS);
//    }
//    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
//    {
//        if (!this->db().objectRegistry::found(rhoName_))
//        {
//            // Rho not available, do not update
//            return;
//        }

//        const fvPatchField<scalar>& rhop =
//            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

//        operator==(*this - n*(n & *this) + n*phip/(rhop*magS));
//    }
//    else
//    {
//        FatalErrorIn
//        (
//            "freeSurfaceVelocityFvPatchVectorField::evaluate()"
//        )
//            << "dimensions of phi are incorrect\n"
//            << "    on patch " << this->patch().name()
//            << " of field " << this->dimensionedInternalField().name()
//            << " in file " << this->dimensionedInternalField().objectPath()
//            << exit(FatalError);
//    }

    fvPatchField<vector>::evaluate();
}


void Foam::freeSurfaceVelocityFvPatchVectorField::write(Ostream& os) const
{
//    fixedGradientFvPatchVectorField::write(os);

    if (phiName_ != "phi")
    {
        os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    }
    if (rhoName_ != "rho")
    {
        os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    }

    os.writeKeyword("secondOrder")
        << secondOrder_ << token::END_STATEMENT << nl;

    os.writeKeyword("nonOrthCorr")
        << nonOrthCorr_ << token::END_STATEMENT << nl;

//    writeEntry("value", os);
    fixedGradientFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        freeSurfaceVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
