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

#include "freeSurfacePressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "fvCFD.H"
#include "faCFD.H"
#include "faMesh.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::freeSurfacePressureFvPatchScalarField::
freeSurfacePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
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
    mu_
    (
        transportProperties_.lookup("mu")
    ),
    rho_
    (
        transportProperties_.lookup("rho")
    ),
    curTimeIndex_(-1)
{}


Foam::freeSurfacePressureFvPatchScalarField::
freeSurfacePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
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
    mu_
    (
        transportProperties_.lookup("mu")
    ),
    rho_
    (
        transportProperties_.lookup("rho")
    ),
    curTimeIndex_(-1)
{
    fvPatchScalarField::operator=(patchInternalField());

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


Foam::freeSurfacePressureFvPatchScalarField::
freeSurfacePressureFvPatchScalarField
(
    const freeSurfacePressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    nonOrthCorr_(ptf.nonOrthCorr_),
    secondOrder_(ptf.secondOrder_),
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
    sigma_(ptf.sigma_),
    sigmaPtr_(ptf.sigmaPtr_),
    mu_(ptf.mu_),
    rho_(ptf.rho_),
    curTimeIndex_(-1)
{}


Foam::freeSurfacePressureFvPatchScalarField::
freeSurfacePressureFvPatchScalarField
(
    const freeSurfacePressureFvPatchScalarField& fcvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fcvpvf, iF),
    nonOrthCorr_(fcvpvf.nonOrthCorr_),
    secondOrder_(fcvpvf.secondOrder_),
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
    sigma_(fcvpvf.sigma_),
    sigmaPtr_(fcvpvf.sigmaPtr_),
    mu_(fcvpvf.mu_),
    rho_(fcvpvf.rho_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::freeSurfacePressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

////    if (curTimeIndex_ != this->db().time().timeIndex())
//    {
//        curTimeIndex_ = this->db().time().timeIndex();

//        label patchID = this->patch().index();

//        const fvMesh& mesh =
//            this->patch().boundaryMesh().mesh();

//        const volVectorField& U =
//            mesh.objectRegistry::lookupObject<volVectorField>("U");

//        // Create Finite Area mesh
//        faMesh aMesh(mesh);

//        // Get surface areas
//        const scalarField& S = aMesh.S();

//        // Evaluate interface normal and get curvature
//        vectorField nf = this->patch().nf();

//        vectorField n = 
//            -S*fac::edgeIntegrate
//            (
//                aMesh.Le()*aMesh.edgeLengthCorrection()
//            )().internalField();

//        forAll (n, faceI)
//        {
//            if (mag(n[faceI])>SMALL)
//            {
//                n[faceI] /= mag(n[faceI]);
//            }
//            else
//            {
//                n[faceI] = nf[faceI];
//            }
//        }

//        // Compute normal component of patch-normal velocity gradient
//        tmp<volVectorField> tUnFs
//        (
//            new volVectorField
//            (
//                IOobject
//                (
//                    "UnFs",
//                    db().time().timeName(),
//                    mesh,
//                    IOobject::NO_READ,
//                    IOobject::NO_WRITE
//                ),
//                mesh,
//                vector::zero
//            )
//        );
//        volVectorField& UnFs = tUnFs();

//        UnFs.boundaryField()[patchID] = n*(n & U.boundaryField()[patchID]);

//        UnFs.boundaryField()[patchID].patchInternalField() = 
//            n*(n & U.boundaryField()[patchID].patchInternalField());

//        string schemeSpec("skewCorrected");
//        tmp<fv::snGradScheme<vector> > tschemeV(
//            fv::snGradScheme<vector>::New(
//                mesh,
//                IStringStream(schemeSpec)()
//            )
//        );
//        fv::snGradScheme<vector>& schemeV = tschemeV();

//        scalarField nGradUn = 
//        (
//            n & schemeV.snGrad(UnFs)().boundaryField()[patchID]
//        );

//        // Evaluate viscous and surface tension forces
////        IOdictionary transportProperties
////        (
////            IOobject
////            (
////                "surfaceProperties",
////                mesh.time().constant(),
////                mesh,
////                IOobject::MUST_READ,
////                IOobject::NO_WRITE
////            )
////        );
////        dimensionedScalar muFluidA(transportProperties.lookup("muFluidA"));
////        dimensionedScalar rhoFluidA(transportProperties.lookup("rhoFluidA"));
////        dimensionedScalar nu = muFluidA/rhoFluidA;


////        const areaScalarField surfaceTension = 
////            mesh.lookupObject<areaScalarField>("surfaceTension");

//        dimensionedScalar muFluidA(transportProperties_.lookup("mu"));
//        dimensionedScalar rhoFluidA(transportProperties_.lookup("rho"));
//        dimensionedScalar nu = muFluidA/rhoFluidA;

//        if 
//        (
//            mesh.foundObject<areaScalarField>("surfaceTension")
//        )
//        {
//            // contaminated interface
//            sigmaPtr_ = const_cast<areaScalarField*> 
//                (&mesh.lookupObject<areaScalarField>("surfaceTension"));
//        }
//        else
//        {
//            // clean interface
//            sigmaPtr_ = new areaScalarField
//            (
//                IOobject
//                (
//                    "sigma",
//                    this->db().time().timeName(),
//                    aMesh.thisDb(),
//                    IOobject::NO_READ,
//                    IOobject::NO_WRITE
//                ),
//                aMesh,
//                sigma_,
//                zeroGradientFaPatchVectorField::typeName
//            );
//        }

//        areaScalarField& surfaceTension = *sigmaPtr_;

//        vectorField surfTensionForce(n.size(), vector::zero);

//        surfTensionForce =
//            S*fac::edgeIntegrate
//            (
//                linearEdgeInterpolate(surfaceTension)
//                *aMesh.Le()*aMesh.edgeLengthCorrection()
//            )().internalField();

//        scalarField normalSurfaceTensionForceDensity = 
//        (
//            n & surfTensionForce
//        )/S;

//        // Set patch-value
//        fvPatchScalarField::operator==
//        (
//            2.*muFluidA.value()*nGradUn
//          - normalSurfaceTensionForceDensity
//        );
//    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    faMesh aMesh(mesh);

    const areaVectorField& nf = aMesh.faceAreaNormals();

    // Get interfacial curvature
    const areaScalarField& K = aMesh.faceCurvatures();

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

    areaScalarField divUs
    (
        fac::div(Us)// - K*(nf&Us) + aMesh.faceCurvatures()*(nf&Us)
    );
//    divUs.correctBoundaryConditions();

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

    // gravity term
    vector pRefPoint(mesh.solutionDict().subDict("PIMPLE").lookup("pRefPoint"));

    // Set patch-value
    fvPatchScalarField::operator==
    (
        -2.0*mu_.value()*divUs.internalField()
      - sigma.internalField()*K.internalField()
//      + rho_.value()
//        *(
//            (
//                mesh.C().boundaryField()[this->patch().index()] 
//              - pRefPoint
//            ) & g_.value()
//        )
   );

    fixedValueFvPatchScalarField::updateCoeffs();
}


Foam::tmp<Foam::Field<Foam::scalar> > 
Foam::freeSurfacePressureFvPatchScalarField::snGrad() const
{
    word pName = this->dimensionedInternalField().name();

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = delta - n*(n&delta);

    vectorField gradpp(this->patch().size(), vector::zero);

    if 
    (
        this->patch().boundaryMesh().mesh()
        .foundObject<volVectorField>("grad(" + pName + ")")
    )
    {
        const fvPatchField<vector>& gradp =
            patch().lookupPatchField<volVectorField, vector>
            (
                "grad(" + pName + ")"
            );

        gradpp = gradp.patchInternalField();
    }
    else
    {
        const volScalarField& p = this->db().lookupObject<volScalarField>(pName);

//        volVectorField gradp("grad(p)", fvc::grad(p));

//        gradpp = patch().lookupPatchField<volVectorField, vector>("grad(p)");

        gradpp = patch().patchField<volVectorField, vector>(fvc::grad(p));

//        gradpp = fvc::grad(p)().boundaryField()[this->patch().index()];
//        vectorField gradpp(this->patch().size(), vector::zero);
    }

    scalarField dpP(this->patch().size(), 0);
    if (nonOrthCorr_)
    {
        dpP = (k & gradpp);
    }

    if (secondOrder_)
    {
        scalarField nGradpP = (n&gradpp);

        return 
            2
           *(
                *this 
              - (this->patchInternalField() + dpP)
            )*this->patch().deltaCoeffs()
          - nGradpP;
    }

    return 
    (
        *this 
      - (patchInternalField() + dpP)
    )*this->patch().deltaCoeffs();
}


Foam::tmp<Foam::Field<Foam::scalar> > 
Foam::freeSurfacePressureFvPatchScalarField::gradientBoundaryCoeffs() const
{
    word pName = this->dimensionedInternalField().name();

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = delta - n*(n&delta);

    vectorField gradpp(this->patch().size(), vector::zero);

    if 
    (
        this->patch().boundaryMesh().mesh()
        .foundObject<volVectorField>("grad(" + pName + ")")
    )
    {
        const fvPatchField<vector>& gradp =
            patch().lookupPatchField<volVectorField, vector>
            (
                "grad(" + pName + ")"
            );

        gradpp = gradp.patchInternalField();
    }
    else
    {
        const volScalarField& p = this->db().lookupObject<volScalarField>(pName);

        gradpp = patch().patchField<volVectorField, vector>(fvc::grad(p));
    }

    scalarField dpP(this->patch().size(), 0);
    if (nonOrthCorr_)
    {
        dpP = (k&gradpp);
    }

    if (secondOrder_)
    {
        scalarField nGradpP = (n&gradpp);

        return
            this->patch().deltaCoeffs()
           *(
                2*(*this - dpP) 
              - this->patchInternalField()
            )
          - nGradpP;
    }

    return this->patch().deltaCoeffs()*(*this - dpP);
}


void Foam::freeSurfacePressureFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
//     writeEntry("value", os);

    os.writeKeyword("secondOrder")
        << secondOrder_ << token::END_STATEMENT << nl;

    os.writeKeyword("nonOrthCorr")
        << nonOrthCorr_ << token::END_STATEMENT << nl;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        freeSurfacePressureFvPatchScalarField
    );
}

// ************************************************************************* //
