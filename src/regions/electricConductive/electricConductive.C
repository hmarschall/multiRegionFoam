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
#include "electricConductive.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

#include "sigmaModelList.H"
#include "dissolvedModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(electricConductive, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        electricConductive,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::electricConductive::electricConductive
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
    phi_
    (
        IOobject
        (
            "phi",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh()
    ),
    j_
    (
        IOobject
        (
            "J",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("J", dimensionSet(0, -3, 0, 0, 0, 1, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    sigmaField_
    (
        IOobject
        (
            "sigma",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("sigma", dimensionSet(-1, -3, 3, 0, 0, 2, 0), 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    T_
    (
        IOobject
        (
            "T",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("T", dimTemperature, 353.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    relax_(electricProperties_.lookupOrDefault<scalar>("relax", 0.0)),
    active_(electricProperties_.lookupOrDefault<Switch>("active", false)),
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
    )

{

    const dictionary& dissolvedDict = cellProperties_.subDict("dissolved");

    if (dissolveOnOff_)
    {
        dissolved_ = dissolvedModel::New(this->mesh(), dissolvedDict);
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

    if (active_)
    {
        active_ = true;
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

Foam::regionTypes::electricConductive::~electricConductive()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::electricConductive::correct()
{
    sigma_->correct(sigmaField_);
    sigmaField_.correctBoundaryConditions();

    if (phi_.needReference() || active_)
    {
        //- Update the potential field.
        //- ElectroNatural, the total electron+proton flux should be zero
        const scalarField& source = j_;
        const scalarField& volume = this->mesh().V();
        scalarField sum = source*volume;
        scalar iDot(Foam::gSum(sum));

        if (active_)
        {
            label patchID = this->mesh().boundaryMesh().findPatchID(patchName_);

            if (patchID == -1)
            {
                FatalErrorInFunction
                    << "Cannot find " << patchName_
                    << "please check the patchName" << exit(FatalError);
            }

            scalarField phiBoundary = phi_.boundaryField()[patchID];

            scalar ibar0 = iDot/Foam::gSum(this->mesh().magSf().boundaryField()[patchID]);

            phi_.boundaryField()[patchID] == phiBoundary + relax_*(ibar0 - ibar_);

            Info << "ibar: " << ibar0 << endl;
            Info << "voltage: " << phiBoundary[0] << endl;
        }
        else
        {
            scalarField& phi = phi_;
            phi -= iDot*relax_;
        }
    }

    if (dissolved_.valid())
    {
        dissolved_->correct();
    }
}

Foam::scalar Foam::regionTypes::electricConductive::getMinDeltaT()
{
    return GREAT;
}

void Foam::regionTypes::electricConductive::setCoupledEqns()
{
}

void Foam::regionTypes::electricConductive::solveRegion()
{
    tmp<fvScalarMatrix> phiEqn
    (
      - fvm::laplacian(sigmaField_, phi_)
      - j_
    );

    phiEqn->relax();

    phiEqn->solve();

    i_ = -sigmaField_ * fvc::grad(phi_);
    i_.correctBoundaryConditions();

    if (dissolved_.valid())
    {
        dissolved_->solve();
    }
}

void Foam::regionTypes::electricConductive::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::electricConductive::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::electricConductive::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::electricConductive::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::electricConductive::meshMotionCorrector()
{
    // do nothing, add as required
}


// ************************************************************************* //
