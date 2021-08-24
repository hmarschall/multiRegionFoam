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

#include "simpleControl.H"


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
		//*this,
        //dimensionedScalar("phi", dimVelocity*dimArea, 0.)
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
//            this->time().timeName(),
//            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    k_(transportProperties_.lookup("k")),
    cp_(transportProperties_.lookup("cp")),
    rhoFluid_(transportProperties_.lookup("rhoFluid")),
    muFluid_(transportProperties_.lookup("muFluid")),
    T_
    (
        IOobject
        (
            "T",
            this->time().timeName(),
            *this,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        *this
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
    )
    
// left from createFields    
//#   include "createUf.H"
//#   include "createSf.H"
//#   include "setRefCell.H"
//#   include "setFluxRequired.H" 
     
/*    
for (label i = 0; i<p.boundaryField().size(); i++)
    {
        if (p.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    },
*/  
//gradp.checkIn(), 
//gradU.checkIn(),     

{}  

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


void Foam::regionTypes::navierStokesFluid::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::navierStokesFluid::solveCoupledPartitioned()
{

/*
    // --- SIMPLEC loop
    label oCorr = 0;
    do
    {
        // Make the fluxes relative
        phi -= fvc::meshPhi(rho, U);

#       include "CourantNo.H"

        //U.storePrevIter();
    
        // Convection-diffusion matrix
        fvVectorMatrix HUEqn
        (
            fvm::div(fvc::interpolate(rho_)*phi_, U_, "div(phi_,U_)")
          - fvm::laplacian(mu_, U_)
        );

        // Time derivative matrix
        fvVectorMatrix ddtUEqn(fvm::ddt(rho_, U_)); 
        
        // solve momentum equation
        solve(ddtUEqn + HUEqn == -gradp_);  
        
        // --- PISO loop
            label corr = 0;
            do
            {
                //p.storePrevIter();
                
                
                
                //- Non-orthogonal pressure corrector loop 
                
                
                
#               include "movingMeshContinuityErrs.H"                

                corr++;
            } while (innerResidual > innerTolerance && corr < nCorr); 
            
            
#           include "surfaceContinuityErrs.H"

#           include "updateSf.H"
#           include "updateUf.H"

            // Check convergence
#           include "checkPIMPLEResidualConvergence.H"

            oCorr++;

        } while (!outerLoopConverged); //- end SIMPLEC loop 
*/                          
}



/*
{
    Info << nl << "Solving ... in " << regionName_ << endl;
    simpleControl simpleControlRegion(*this);

    dimensionedScalar alpha = k_/(rho_*cp_);

    while (simpleControlRegion.correctNonOrthogonal())
    {
        tmp<fvScalarMatrix> TEqn
        (
            fvm::ddt(T_)
          + fvm::div(phi_, T_)
         ==
            fvm::laplacian(alpha, T_)
        );

        TEqn->relax();

        TEqn->solve();
    }
}
*/



// ************************************************************************* //
