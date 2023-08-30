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

#include "surfaceStandardDWModel.H"
#include "regionType.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceDissolvedModels
{
    defineTypeNameAndDebug(surfaceStandardDWModel, 0);

    addToRunTimeSelectionTable
    (
        surfaceDissolvedModel,
        surfaceStandardDWModel,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Private functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::surfaceDissolvedModels::surfaceStandardDWModel::DNaf
(
    const scalar& d0,
    const scalar& w
) const
{
    scalar tmp(0.0);

    if(w <= 3.0)
    {
        tmp = d0*3.1e-7*w*(pow(M_E, 0.28*w) - 1.0);
    }
    else
    {
        tmp = d0*4.17e-8*w*(161*pow(M_E, -w) + 1.0);
    }

    return tmp;
}


Foam::scalar Foam::surfaceDissolvedModels::surfaceStandardDWModel::effLam(const scalar& act) const
{
    scalar tmp(0.0);

    if(act <= 1.0)
    {
        tmp = 0.043 + 17.81*act - 39.85*Foam::pow(act, 2) + 36.0*Foam::pow(act, 3);
    }
    else if (act > 3.0)
    {
        tmp = 15.8;
    }
    else
    {
        tmp = 14.0 + 1.4*(act - 1.0);
    }

    return tmp;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceDissolvedModels::surfaceStandardDWModel::surfaceStandardDWModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    surfaceDissolvedModel(mesh),

    dict_(dict),
    ksi_(dict_.lookup("ksi")),
    nd_("nd", dimless, readScalar(dict_.lookup("nd"))),
    rhoOnEW_("rhoOnEW", dimMoles/dimVol, readScalar(dict_.lookup("rhoOnEW"))),
    iName_(dict_.lookupOrDefault<word>("i", "i")),
    TName_(dict_.lookupOrDefault<word>("TElectric", "TElectric")),
    relax_(dict_.lookupOrDefault<scalar>("relax", 0.0)),
    corr_(dict_.lookupOrDefault<scalar>("corr", 1.0))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceDissolvedModels::surfaceStandardDWModel::~surfaceStandardDWModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::surfaceDissolvedModels::surfaceStandardDWModel::solve()
{
    const volVectorField& i = mesh().lookupObject<volVectorField>(iName_);

    surfaceScalarField phiI = fvc::interpolate(i) & mesh().Sf();

    tmp<fvScalarMatrix> wEqn
    (
        fvm::ddt(rhoOnEW_, lambda_)
      + fvm::div(phiI*nd_/Foam::phaseModel::dimF, lambda_, "div(D,lambda)")
      - fvm::laplacian(rhoOnEW_*Dwm_, lambda_, "laplacian(diff,lambda)")
      - dmdt_
    );

    //- Set reference values
    {
        const scalarField& source = dmdt_;
        const scalarField& volume = mesh().V();

        // Update the water content
        // Water intake = water uptake
        scalarField sum = source*volume;
        scalar iDot(Foam::gSum(sum)/Foam::gSum(volume));

        for (label i = 0; i < lambda_.size(); i++)
        {
            wEqn->setReference(i, lambda_[i] + iDot*relax_);
        }
    }

    wEqn->relax();
    wEqn->solve();
}


void Foam::surfaceDissolvedModels::surfaceStandardDWModel::correct()
{
    scalarField& lambdaIn = lambda_;

    const scalarField& T = mesh().lookupObject<volScalarField>(TName_);

    scalarField Dwm0 = Foam::exp(-2436/T)*corr_;

    forAll(Dwm_, cellI)
    {
        Dwm_[cellI] = DNaf(Dwm0[cellI], lambdaIn[cellI]);
    }

    Dwm_.correctBoundaryConditions();

    // The activity on the boundary is updated in electrochemicalReaction
    // for the interface case

    // For the interface case:
    // I will map the calculated act_ (from electrochemicalReactions)
    // onto the face cells in this region
    // >> Dann kann ich ja ruhig correctBC hier machen, da es ja eh eine zeroGradient BC ist
    // if(!this->surfaceElectrochemistry_)
    // {
        // Allerdings brauche ich act nur für den update des übertragenen Wassers
        // d.h. ich brauche es eig. nur auf der Boundary um den Wert davon abzugreifen
        //  act_.correctBoundaryConditions();
    // }

    const scalarField& source = dmdt_;
    const scalarField& volume = mesh().V();

    // Update the water content
    // Water intake = water uptake
    scalarField sum = source*volume;
    scalar iDot(Foam::gSum(sum)/Foam::gSum(volume));
    this->lambda_.internalField() += iDot*relax_;
}


void Foam::surfaceDissolvedModels::surfaceStandardDWModel::update
(
    const word& clName
)
{

    label ionCatalystID = mesh().boundaryMesh().findPatchID(clName);
    const polyPatch& ionCatalystPatch = mesh().boundaryMesh()[ionCatalystID];

    const fvPatchField<scalar>& lambdaInPatch = lambda_.boundaryField()[ionCatalystID];
    const volScalarField& lambdaIn = lambda_;
    fvPatchField<scalar>& dmdtInPatch = dmdt_.boundaryField()[ionCatalystID];
    volScalarField& dmdtIn = dmdt_;
    const fvPatchField<scalar>& actInPatch = this->act_.boundaryField()[ionCatalystID];
    scalarField actInPatchSf = actInPatch;
    // Info << "Print actInPatchSf in dissolvedModel " << actInPatchSf << endl;
    const volScalarField& actIn = this->act_;
    volScalarField& lambda = lambda_;
    const fvPatchField<scalar>& DwmIn = Dwm_.boundaryField()[ionCatalystID];

    // Gleichung passt aber weder zu den Einheiten im Volumen [mol/m^3/s] noch zu [mol/m^2/s]
    // Da diese Gleichung allerdings davor für die Berechnung von dmdt im Volumen benutzt wurde,
    // muss man für die cell values mal V geteilt durch A machen

    // Kann ich ja hier hin machen und dann die faces updaten
    dmdt_.correctBoundaryConditions();

    // Only consider face cells for water adsorption or desorption
    // Smear to face cells
    forAll(ionCatalystPatch, faceI)
    {
        // For act only the boundary is filled with values in electrochemicalReactions
        label faceCelli = ionCatalystPatch.faceCells()[faceI];
        dmdtIn[faceCelli] =
            ksi_[clName]*rhoOnEW_.value()
            * (effLam(actInPatch[faceI]) - lambdaIn[faceCelli]);

        dmdtInPatch[faceI] = dmdtIn[faceCelli] * (this->mesh().V()[faceCelli] / dmdtInPatch.patch().magSf()[faceI]); 
    }
}


bool Foam::surfaceDissolvedModels::surfaceStandardDWModel::read(const dictionary& dict)
{
    const dictionary& dict0 = dict.subDict(type() + "Coeffs");

    dict0.lookup("ksi") >> ksi_;
    nd_ = dimensionedScalar("nd", dimless, readScalar(dict0.lookup("nd")));
    rhoOnEW_ = dimensionedScalar("rhoOnEW", dimMoles/dimVol, readScalar(dict0.lookup("rhoOnEW")));
    iName_ = word(dict0.lookupOrDefault<word>("i", "i"));
    TName_ = word(dict0.lookupOrDefault<word>("T", "T"));
    relax_ = scalar(dict0.lookupOrDefault<scalar>("relax", 0.0));
    corr_ = scalar(dict0.lookupOrDefault<scalar>("corr", 1.0));

    return true;
}
// ************************************************************************* //
