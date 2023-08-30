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
#include "surfaceElectricConductive.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

#include "sigmaModelList.H"
#include "surfaceDissolvedModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(surfaceElectricConductive, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        surfaceElectricConductive,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::surfaceElectricConductive::surfaceElectricConductive
(
    const Time& runTime,
    const word& regionName
)
:
    regionType(runTime, regionName),

    regionName_(regionName),

    electricProperties_
    (
        IOobject
        (
            "electricProperties",
            this->mesh().time().constant(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    i_
    (
        IOobject
        (
            "i",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedVector("i", dimensionSet(0, -2, 0, 0, 0, 1, 0), vector::zero),
        zeroGradientFvPatchVectorField::typeName
    ),
    T_
    (
        IOobject
        (
            "TElectric",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("TElectric", dimTemperature, 353.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    relax_(electricProperties_.lookupOrDefault<scalar>("relax", 0.0)),
    control_(electricProperties_.lookupOrDefault<Switch>("active", false)),
    dissolveOnOff_(electricProperties_.lookupOrDefault<Switch>("dissolveOnOff", false)),
    patchName_(word::null),
    galvanostatic_(true),
    ibar_(0.0),
    voltage_(0.0),
    cellProperties_
    (
        IOobject
        (
            "cellProperties",
            this->mesh().time().constant(),        // global mesh
            this->mesh().time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    sigmaField_(nullptr),
    phi_(nullptr)
{

    // I need sigma also updated at the boundary
    // Dunno how to set the zeroGrad BC otherwise for now
    sigmaField_ = lookupOrRead<volScalarField>(mesh(), "sigma");


        // set potential field
    phi_ = lookupOrRead<volScalarField>(mesh(), "phi");

    const dictionary& dissolvedDict = cellProperties_.subDict("dissolved");

    if (dissolveOnOff_)
    {
        dissolved_ = surfaceDissolvedModel::New(this->mesh(), dissolvedDict);
    }

    sigma_.set
    (
        new sigmaModelList
        (
            this->mesh(),
            electricProperties_.subDict("sigma")
        )
    );

    const dictionary& galvanostaticDict = cellProperties_.subDict("galvanostatic");

    if (control_)
    {
        control_ = true;
        galvanostatic_ = Switch(galvanostaticDict.lookupOrDefault("active", true));
        patchName_ = word(galvanostaticDict.lookup("patchName"));

        if (galvanostatic_)
        {
            galvanostaticDict.lookup("ibar") >> ibar_;
        }
        else
        {
            galvanostaticDict.lookup("voltage") >> voltage_;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::surfaceElectricConductive::~surfaceElectricConductive()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::surfaceElectricConductive::correct()
{
    sigma_->correct(sigmaField_());
    sigmaField_->correctBoundaryConditions();

    // Since the phiEqn are coupled, need to update the current density here
    i_ = -sigmaField_() * fvc::grad(phi_());
    i_.correctBoundaryConditions();

    if (phi_().needReference() || control_)
    {
        if (control_)
        {
            label patchID = this->mesh().boundaryMesh().findPatchID(patchName_);
            fvPatchField<scalar>& phiPatch = phi_().boundaryField()[patchID];
            const fvPatchField<scalar>& sigmaPatch = sigmaField_->boundaryField()[patchID];
            const fvPatchField<vector>& iPatch = i_.boundaryField()[patchID];

            if (patchID == -1)
            {
                FatalErrorInFunction
                    << "Cannot find " << patchName_
                    << "please check the patchName" << exit(FatalError);
            }

            // change the fixedValue of phi on interconnect0
            scalarField phiBoundary = phi_().boundaryField()[patchID];

            // // Calculate the average of sigma on this patch (for the case sigma is not constant there)
            const scalar avgSigmaPatch = Foam::gAverage(sigmaPatch);
            scalarField sum = iPatch & iPatch.patch().Sf();
            scalar ibar0 = Foam::gSum(sum)/Foam::gSum(iPatch.patch().magSf());

            phi_().boundaryField()[patchID] == phiBoundary + relax_*(ibar0 - ibar_);

            Info << "ibar: " << ibar0 << "\t"<< "voltage: " << Foam::gAverage(phiBoundary) << endl;
        }
    }

    if (dissolved_.valid())
    {
        dissolved_->correct();
    }
}

Foam::scalar Foam::regionTypes::surfaceElectricConductive::getMinDeltaT()
{
    return GREAT;
}

void Foam::regionTypes::surfaceElectricConductive::setCoupledEqns()
{
    Info << "setCoupledEqns in region " << this->name() << endl;

    phiEqn =
    (
        - fvm::laplacian(sigmaField_(), phi_())
    );

    fvScalarMatrices.set
    (
        phi_().name()
      + mesh().name() + "Mesh"
      + surfaceElectricConductive::typeName + "Type"
      + "Eqn",
        &phiEqn()
    );
}

void Foam::regionTypes::surfaceElectricConductive::solveRegion()
{
    if (dissolved_.valid())
    {
        dissolved_->solve();
    }
}

void Foam::regionTypes::surfaceElectricConductive::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::surfaceElectricConductive::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::surfaceElectricConductive::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::surfaceElectricConductive::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::surfaceElectricConductive::meshMotionCorrector()
{
    // do nothing, add as required
}


// ************************************************************************* //
