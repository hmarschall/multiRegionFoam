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
#include "solutionControl.H"

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
    const Time& runTime,
    const word& regionName
)
:
    regionType(runTime, regionName),
    aMesh_(*this),
    
    regionName_(regionName),

    U_
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
    ),
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
		linearInterpolate(U_) & (*this).Sf()    
    ),
    p_
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
    ),
    
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
    rhoFluid_
    (
        // transportProperties_.subDict(regionName_).lookup("rho")
            //currently there is no subdict for region name, each region has 
            //a transportProperties dictionary in its folder
        transportProperties_.lookup("rho")
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
    muFluid_
    (
        //transportProperties_.subDict(regionName_).lookup("mu")
        transportProperties_.lookup("mu")
    ),
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
		fvc::grad(p_)   
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
		fvc::grad(U_) 
    ), 
    pcorrTypes_
    (
        p_.boundaryField().size(),
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
        dimensionedScalar("pcorr", p_.dimensions(), 0.0),
        pcorrTypes_
    ),
    UUrf_ 
    ( 
        this->solutionDict().subDict("relaxationFactors")
        .lookupOrDefault<scalar>(U_.name(), 1)   
    ), 
    
    pRefCell_(0),
    pRefValue_(0),     
    
    innerResidual_(1),
    residualPressure_(1),

    corr_(0),
    corrNonOrtho_(0)
{
    
    for (label i = 0; i<p_.boundaryField().size(); i++)
    {
        if (p_.boundaryField()[i].fixesValue())
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
    // do nothing, add as required
}


void Foam::regionTypes::navierStokesFluid::updateFields()
{
    // do nothing, add as required
}


void Foam::regionTypes::navierStokesFluid::solveRegion()
{
    // --- PIMPLE loop      
    pimpleControl pimple(*this);
         
    while (pimple.loop())
    // (corr_ == nCorrPIMPLE_ + 1) 
    {
        U_.storePrevIter();     
        
        // Convection-diffusion matrix
        fvVectorMatrix HUEqn
        (
            fvm::div(fvc::interpolate(rho_)*phi_, U_, "div(phi,U)")
          - fvm::laplacian(mu_, U_)
        );

        // Time derivative matrix
        fvVectorMatrix ddtUEqn(fvm::ddt(rho_, U_));

        // UEqn
        fvVectorMatrix UEqn(ddtUEqn + HUEqn);
        
        // Save source and boundaryCoeffs
        vectorField S0 = UEqn.source();
        FieldField<Field, vector> B0 = UEqn.boundaryCoeffs();
        
        // Relax and solve momentum equation
        UEqn.relax(UUrf_);

        if (pimple.momentumPredictor()) 
        {    
            solve(UEqn == -gradp_);
        }
        
        // reset equation to ensure relaxation parameter 
        // is not causing problems with time consistency
        UEqn = (ddtUEqn + HUEqn);
        UEqn.source() = S0;
        UEqn.boundaryCoeffs() = B0;

        // --- PISO loop
        while (pimple.correct())  
        // (corrPISO_ <= nCorrPISO_)            
        {
            p_.storePrevIter();   
                 
            AU_ = UEqn.A();
            HU_ = UEqn.H();

            volVectorField HbyA("HbyA", U_);

            HbyA = HU_/AU_;
            
            surfaceScalarField phiHbyA // moved from the constructor to avoid 
                                       // error with divide operator at runtime
            (
                "phiHbyA",
                (
                    (fvc::interpolate(HU_)/fvc::interpolate(AU_))
                    & this->Sf()
                )
            );

            phiHbyA += fvc::ddtPhiCorr((1.0/AU_)(), rho_, U_, phiHbyA);

            tmp<volScalarField> AtU(AU_);
            AtU = max(AU_ - UEqn.H1(), 0.1*AU_);

            surfaceScalarField AtUf("AtUf", fvc::interpolate(AtU()));
            surfaceScalarField AUf("AUf", fvc::interpolate(AU_));

            phiHbyA += (1.0/AtUf - 1.0/AUf)*fvc::snGrad(p_)*this->magSf();

            HbyA -= (1.0/AU_ - 1.0/AtU)*gradp_;

            forAll(phiHbyA.boundaryField(), patchI)
            {
                if
                (
                    !p_.boundaryField()[patchI].fixesValue()
                 && isA<zeroGradientFvPatchVectorField>
                    (
                        U_.boundaryField()[patchI]
                    )
                )
                {
                    phiHbyA.boundaryField()[patchI] =
                    (
                        U_.boundaryField()[patchI]
                        & this->Sf().boundaryField()[patchI]
                    );
                }
            }      
                 
            while (pimple.correctNonOrthogonal())
            // (corrNonOrtho_ <= nNonOrthCorr_ + 1)
            // vs. for (label nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(1.0/AtUf, p_) == fvc::div(phiHbyA)
                );
                
                // #include "setRefCell.H" 
                    // error: ‘interface’ was not declared
                                          
                // #include "setReference.H" // need "setRefCell.H"   
                
                pEqn.setReference(pRefCell_, pRefValue_);             
                
                pEqn.solve
                (
                    this->solutionDict().solver
                    (
                        p_.select(pimple.finalInnerIter())
                    )
                ); 
                 

                innerResidual_ = pEqn.solve().initialResidual();

                if (corrNonOrtho_ == 0 && corr_ == 0) 
                {
                    residualPressure_ = innerResidual_;
                }  
                        
                if (pimple.finalNonOrthogonalIter()) 
                // (corrNonOrtho_ == nNonOrthCorr_ + 1)
                {
                    phi_ = phiHbyA - pEqn.flux();
                }                            
            }
            
            //- Pressure relaxation except for last corrector
            if (!pimple.finalIter())
            {
                p_.relax();
            }

            // Update of pressure gradient
            gradp_ = fvc::grad(p_);

            // Momentum corrector
            U_ = UUrf_*(HbyA - 1.0/AtU*gradp_ + (1 - UUrf_)*U_.prevIter());

            U_.correctBoundaryConditions();

            // Update of velocity gradient
            gradU_ = fvc::grad(U_);
                
        }
        //} while (innerResidual > innerTolerance && corr < nCorr); 
    }
}

// ************************************************************************* //
