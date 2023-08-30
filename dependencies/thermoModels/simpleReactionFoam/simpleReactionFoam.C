/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

Application
    rhoReactingFoam

Description
    Solver for combustion with chemical reactions using density based
    thermodynamics package.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoReactionThermo.H"
#include "psiReactionThermo.H"
#include "rhoThermo.H"
#include "multiComponentMixture.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    pimpleControl pimple(mesh);

// #   include "createFields.H"
Info<< nl << "Reading thermophysicalProperties" << endl;

autoPtr<rhoReactionThermo> thermo = rhoReactionThermo::New(mesh);

Info<< nl << "Print 1" << endl;

basicSpecieMixture& composition = thermo->composition();

Info<< nl << "Print 2" << endl;

PtrList<volScalarField>& Y = thermo->composition().Y();

forAll(Y, i)
{
    Info << "Print Y[i].name() " << Y[i].name() << endl;
    Info << "Print Foam::gAverage(Y[i]) " << Foam::gAverage(Y[i]) << endl;
}

// autoPtr<basicThermo> thermo2 = thermo;

// PtrList<volScalarField>& Y2 = thermo2->composition().Y();

// forAll(Y2, i)
// {
//     Info << "Print Y2[i].name() " << Y2[i].name() << endl;
//     Info << "Print Foam::gAverage(Y2[i]) " << Foam::gAverage(Y2[i]) << endl;
// }

// const PtrList<rhoReactionThermo> speciesThermo_(1);
// speciesThermo_ = 
//         dynamic_cast<const multiComponentMixture<rhoReactionThermo>&>
//             (thermo)->speciesData();



// word inertSpecie(thermo->lookup("inertSpecie"));

// volScalarField rho
// (
//     IOobject
//     (
//         "rho",
//         runTime.timeName(),
//         mesh
//     ),
//     thermo->rho()
// );


Info<< nl << "Ending simpleREactionFoam " << endl;
// const volScalarField& psi = thermo->psi();
// volScalarField& h = thermo->h();
// const volScalarField& T = thermo->T();

    
    // For now just get access to the data I want to in the thermo model
    

    return 0;
}


// ************************************************************************* //
