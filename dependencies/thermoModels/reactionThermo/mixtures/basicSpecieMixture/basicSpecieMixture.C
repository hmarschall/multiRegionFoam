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

\*---------------------------------------------------------------------------*/

#include "basicSpecieMixture.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicSpecieMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicSpecieMixture::basicSpecieMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicMixture(thermoDict, mesh, phaseName),
    // phaseName_(phaseName),
    species_(specieNames),
    defaultSpecie_
    (
        species_.size()
      ? word(thermoDict.lookup
        (
            {"defaultSpecie", "inertSpecie"}
        ))
      : word("undefined")
    ),
    defaultSpecieIndex_(-1),
    active_(species_.size(), true),
    Y_(species_.size())
    // X_(species_.size()),
    // C_(species_.size())
{

    // if (species_.size())
    // {
    //     if (species_.contains(defaultSpecie_))
    //     {
    //         defaultSpecieIndex_ = species_[defaultSpecie_];
    //     }
    //     else
    //     {
    //         FatalIOErrorInFunction(thermoDict)
    //             << "default specie " << defaultSpecie_
    //             << " not found in available species " << species_
    //             << exit(FatalIOError);
    //     }
    // }

    forAll(species_, i)
    {
        IOobject header
        (
            // IOobject::groupName("Y" + species_[i], phaseName),
            IOobject::groupName(species_[i], phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        // check if field exists and can be read
        if (header.headerOk())
        {
            Info << "Create vsf Yi for specie " << species_[i] << endl;
            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        // IOobject::groupName("Y" + species_[i], phaseName),
                        IOobject::groupName(species_[i], phaseName),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
            // X_.set
            // (
            //     i,
            //     new volScalarField
            //     (
            //         IOobject
            //         (
            //             IOobject::groupName("X_" + species_[i], phaseName),
            //             mesh.time().timeName(),
            //             mesh,
            //             IOobject::No_READ,
            //             IOobject::AUTO_WRITE
            //         ),
            //         mesh
            //     )
            // );
            // C_.set
            // (
            //     i,
            //     new volScalarField
            //     (
            //         IOobject
            //         (
            //             IOobject::groupName("C_" + species_[i], phaseName),
            //             mesh.time().timeName(),
            //             mesh,
            //             IOobject::NO_READ,
            //             IOobject::AUTO_WRITE
            //         ),
            //         mesh
            //     )
            // );

            Info << "End vsf Yi for specie " << species_[i] << endl;
        }
        else
        {
            volScalarField Ydefault
            (
                IOobject
                (
                    IOobject::groupName("Ydefault", phaseName),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        // IOobject::groupName("Y" + species_[i], phaseName),
                        IOobject::groupName(species_[i], phaseName),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    Ydefault
                )
            );
            // X_.set
            // (
            //     i,
            //     new volScalarField
            //     (
            //         IOobject
            //         (
            //             IOobject::groupName("X_" + species_[i], phaseName),
            //             mesh.time().timeName(),
            //             mesh,
            //             IOobject::NO_READ,
            //             IOobject::AUTO_WRITE
            //         ),
            //         Ydefault
            //     )
            // );
            // C_.set
            // (
            //     i,
            //     new volScalarField
            //     (
            //         IOobject
            //         (
            //             IOobject::groupName("C_" +species_[i], phaseName),
            //             mesh.time().timeName(),
            //             mesh,
            //             IOobject::NO_READ,
            //             IOobject::AUTO_WRITE
            //         ),
            //         Ydefault
            //     )
            // );
        }
    }

    // Do not enforce constraint of sum of mass fractions to equal 1 here
    // - not applicable to all models
}


// ************************************************************************* //
