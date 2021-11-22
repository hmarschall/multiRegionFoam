/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
        
#include "fvCFD.H"
#include "faCFD.H"
#include "navierStokesFluid.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(navierStokesFluid, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        navierStokesFluid,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::navierStokesFluid::navierStokesFluid
(
    const fvMesh& mesh,
    const word& regionName
)
:
    regionType(mesh, regionName),
    mesh_(mesh),
    aMesh_(*this),
    
    regionName_(regionName),

    //! coupled fields changed to ptr
    U_(nullptr), 
    p_(nullptr),  
    
//    U_
//    (
//        IOobject
//        (
//            "U",
//            this->time().timeName(),
//            *this,
//            IOobject::MUST_READ,
//            IOobject::AUTO_WRITE
//        ),
//        *this
//    ),
//    p_
//    (
//		IOobject
//		(
//			"p",
//            this->time().timeName(),
//            *this,
//			IOobject::MUST_READ,
//			IOobject::AUTO_WRITE
//		),
//		*this   
//    ),
    
    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            this->time().constant(),
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    rhoFluid_(transportProperties_.lookup("rhoFluid")),
    muFluid_(transportProperties_.lookup("muFluid")),
    mu_
    (
        IOobject
        (
            "mu",
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        muFluid_
    ),
    rho_
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        rhoFluid_
    ),
    
    surfaceProperties_
    (
        IOobject
        (
            "surfaceProperties",
            this->time().constant(),
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    cleanSurfaceTension_(surfaceProperties_.lookup("cleanSurfaceTension")),
    surfaceTension_
    (
        IOobject
        (
            "surfaceTension",
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_,   
        cleanSurfaceTension_
    ),
    
    AU_
    (
        IOobject
        (
            "AU",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ),
        *this,
        dimMass/(dimVolume*dimTime),
        zeroGradientFvPatchScalarField::typeName
    ),
    HU_
    (
        IOobject
        (
            "HU",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this,
        dimMass/sqr(dimLength*dimTime),
        zeroGradientFvPatchVectorField::typeName
    ),
    phiHbyA_
    (
        "phiHbyA",
        (
            (fvc::interpolate(HU_)/fvc::interpolate(AU_))
            & this->Sf()
        )
    ),

// Fields below depend on U_ or p_ which are nullptr and not set yet -----------
// this yields nullptr error when running a case (change these to ptr also?)
    phi_
    (
		IOobject
		(
			"phi",
            this->time().timeName(),
            *this,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		linearInterpolate(U_()) & (*this).Sf()    
    ),
    
    UUrf_ // from UEqn.H 
    ( 
        this->solutionDict().subDict("relaxationFactors")
        .lookupOrDefault<scalar>(U_().name(), 1)   
    ), 
    
    gradp_
    (
		IOobject
		(
			"grad(p)",
            this->time().timeName(),
            *this,
			IOobject::NO_READ,
            IOobject::NO_WRITE
		),
		fvc::grad(p_())   
    ), 
    gradU_
    (
		IOobject
		(
			"grad(U)",
            this->time().timeName(),
            *this,
			IOobject::NO_READ,
            IOobject::NO_WRITE
		),
		fvc::grad(U_()) 
    ), 
    pcorrTypes_
    (
        p_().boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    ),   
    pcorr_
    (
        IOobject
        (
            "pcorr",
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        dimensionedScalar("pcorr", p_().dimensions(), 0.0),
        pcorrTypes_
    )
{

    U_.reset
    (
        new volVectorField
        (
           IOobject
           (
               "U",
               this->time().timeName(),
               *this,
               IOobject::MUST_READ,
               IOobject::AUTO_WRITE
           ),
           *this
        )
    ); 
    
    p_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "p",
                this->time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        )
    );

    // from createFields.H     
    for (label i = 0; i<p_().boundaryField().size(); i++)
    {
        if (p_().boundaryField()[i].fixesValue())
        {
            pcorrTypes_[i] = fixedValueFvPatchScalarField::typeName;
        }
    };  
}  

// left from createFields    
//#   include "createUf.H"
//#   include "createSf.H"
//#   include "setRefCell.H"
//#   include "setFluxRequired.H" 
     
//gradp.checkIn(), 
//gradU.checkIn(),     

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::navierStokesFluid::~navierStokesFluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::navierStokesFluid::correct()
{
    // do nothing, add as required
}


void Foam::regionTypes::navierStokesFluid::setRDeltaT()
{
    // do nothing, add as required
}


void Foam::regionTypes::navierStokesFluid::setCoupledEqns()
{
    pimpleControl pimple(*this);

// from UEqn.H 
    // Momentum equation
    fvVectorMatrix UEqn = 
    (
        fvm::div(fvc::interpolate(rho_)*phi_, U_(), "div(phi,U)")
      - fvm::laplacian(mu_, U_())
      + fvm::ddt(rho_, U_())
      ==
      -gradp_
    );
    
    // Ignored for now:
    // Save source and boundaryCoeffs
    // vectorField S0 = UEqn.source();
    // FieldField<Field, vector> B0 = UEqn.boundaryCoeffs();

    UEqn.relax(UUrf_);
    
    fvVectorMatrices.set
    (
        U_().name() + this->name() + "Eqn", 
        new fvVectorMatrix(UEqn)
    );
       
    // Ignored for now:
    // Reset equation to ensure relaxation parameter 
    // is not causing problems with time consistency
    // UEqn = (ddtUEqn + HUEqn);
    // UEqn.source() = S0;
    // UEqn.boundaryCoeffs() = B0;

    
    // pressure equation

    // check: UEqn and U_() are called here before or after solving 
    // the momentum equation #361
    // check: U.storePrevIter() & p.storePrevIter() in interTrackFoam 
    // also U_().prevIter() used in momentum corrector #463
    
    AU_ = UEqn.A();
    HU_ = UEqn.H();

    volVectorField HbyA("HbyA", U_());

    HbyA = HU_/AU_;

    phiHbyA_ += fvc::ddtPhiCorr((1.0/AU_)(), rho_, U_(), phiHbyA_);

    tmp<volScalarField> AtU(AU_);
    AtU = max(AU_ - UEqn.H1(), 0.1*AU_);

    surfaceScalarField AtUf("AtUf", fvc::interpolate(AtU()));
    surfaceScalarField AUf("AUf", fvc::interpolate(AU_));

    phiHbyA_ += (1.0/AtUf - 1.0/AUf)*fvc::snGrad(p_())*this->magSf();

    HbyA -= (1.0/AU_ - 1.0/AtU)*gradp_;

    forAll(phiHbyA_.boundaryField(), patchI)
    {
        if
        (
            !p_().boundaryField()[patchI].fixesValue()
         && isA<zeroGradientFvPatchVectorField>
            (
                U_().boundaryField()[patchI]
            )
        )
        {
            phiHbyA_.boundaryField()[patchI] =
            (
                U_().boundaryField()[patchI]
                & this->Sf().boundaryField()[patchI]
            );
        }
    }

//    while (pimple.correctNonOrthogonal())
//    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(1.0/AtUf, p_()) == fvc::div(phiHbyA_)
        );
//    }
    // check: setting pRefCell and pRefValue for #439
    // in icoFoam: 
        // label pRefCell = 0;
        // scalar pRefValue = 0.0;  
    
    // in surfaceTrackingFoam:
        // #include "setReference.H" 
        // also needs #include "setRefCell.H" 
        // but this gives error "interface not decalred"

    // pEqn.setReference(pRefCell, pRefValue);

    fvScalarMatrices.set
    (
        p_().name() + this->name() + "Eqn",
        new fvScalarMatrix(pEqn)
    );

//    // Ignored  (should be in multiRegionSystem?)  
//    //    innerResidual = pEqn.solve().initialResidual();

//    //    if ( nonOrth == 0 && corr == 0 )
//    //    {
//    //        residualPressure = innerResidual;
//    //    }   
// 
//// check: the below belong here or updateField() 
//    //- Pressure relaxation
//    p_().relax();

//    //- Update of pressure gradient
//    gradp_ = fvc::grad(p_());

//    //- Momentum corrector
//    U_() = UUrf_*(HbyA - 1.0/AtU*gradp_ + (1 - UUrf_)*U_().prevIter());

//    U_().correctBoundaryConditions();

//    //- Update of velocity gradient
//    gradU_ = fvc::grad(U_());
}


// TODO:
void Foam::regionTypes::navierStokesFluid::updateFields()
{
    phi_ = phiHbyA_ 
           - (fvScalarMatrices[p_().name() + this->name() + "Eqn"])->flux();
}


void Foam::regionTypes::navierStokesFluid::solveRegion()
{
    // do nothing, add as required
}

// ************************************************************************* //
