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

#include "basicThermo.H"
#include "fvMesh.H"
#include "HashTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedEnthalpyFvPatchScalarField.H"
#include "gradientEnthalpyFvPatchScalarField.H"
#include "mixedEnthalpyFvPatchScalarField.H"
// #include "fixedInternalEnergyFvPatchScalarField.H"
// #include "gradientInternalEnergyFvPatchScalarField.H"
// #include "mixedInternalEnergyFvPatchScalarField.H"

/* * * * * * * * * * * * * * Private Static Data * * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(basicThermo, 0);
    defineRunTimeSelectionTable(basicThermo, fvMesh);
}

const Foam::word Foam::basicThermo::dictName("thermophysicalProperties");

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::basicThermo::heBoundaryTypes()
{
    const volScalarField::GeometricBoundaryField& tbf = T_.boundaryField();

    wordList hbt = tbf.types();

    forAll(tbf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = fixedEnthalpyFvPatchScalarField::typeName;
        }
        else if
        (
            isA<zeroGradientFvPatchScalarField>(tbf[patchi])
         || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
        )
        {
            hbt[patchi] = gradientEnthalpyFvPatchScalarField::typeName;
        }
        else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = mixedEnthalpyFvPatchScalarField::typeName;
        }
    }

    return hbt;
}


void Foam::basicThermo::heBoundaryCorrection(volScalarField& h)
{
    // Info << "Print hbf davor " << h.boundaryField() << endl;
    volScalarField::GeometricBoundaryField& hbf = h.boundaryField();

    forAll(hbf, patchi)
    {
        if (isA<gradientEnthalpyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<gradientEnthalpyFvPatchScalarField>(hbf[patchi]).gradient()
                = hbf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedEnthalpyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<mixedEnthalpyFvPatchScalarField>(hbf[patchi]).refGrad()
                = hbf[patchi].fvPatchField::snGrad();
        }
    }
    // Info << "Print hbf danach " << h.boundaryField() << endl;
}


// Foam::wordList Foam::basicThermo::eBoundaryTypes()
// {
//     const volScalarField::GeometricBoundaryField& tbf = T_.boundaryField();

//     wordList ebt = tbf.types();

//     forAll(tbf, patchi)
//     {
//         if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
//         {
//             ebt[patchi] = fixedInternalEnergyFvPatchScalarField::typeName;
//         }
//         else if
//         (
//             isA<zeroGradientFvPatchScalarField>(tbf[patchi])
//          || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
//         )
//         {
//             ebt[patchi] = gradientInternalEnergyFvPatchScalarField::typeName;
//         }
//         else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
//         {
//             ebt[patchi] = mixedInternalEnergyFvPatchScalarField::typeName;
//         }
//     }

//     return ebt;
// }


// void Foam::basicThermo::eBoundaryCorrection(volScalarField& e)
// {
//     volScalarField::GeometricBoundaryField& ebf = e.boundaryField();

//     forAll(ebf, patchi)
//     {
//         if (isA<gradientInternalEnergyFvPatchScalarField>(ebf[patchi]))
//         {
//             refCast<gradientInternalEnergyFvPatchScalarField>(ebf[patchi])
//                 .gradient() = ebf[patchi].fvPatchField::snGrad();
//         }
//         else if (isA<mixedInternalEnergyFvPatchScalarField>(ebf[patchi]))
//         {
//             refCast<mixedInternalEnergyFvPatchScalarField>(ebf[patchi])
//                 .refGrad() = ebf[patchi].fvPatchField::snGrad();
//         }
//     }
// }

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicThermo::basicThermo(const fvMesh& mesh, const word& phaseName) // obj is the mesh (see hREactionThermo.C)
:
    IOdictionary
    (
        IOobject
        (
            phasePropertyName(dictName, phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    phaseName_(phaseName),
    p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    T_
    (
        IOobject
        (
            phasePropertyName("T"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    alpha_
    (
        IOobject
        (
            phasePropertyName("thermo:alpha"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("alpha", dimensionSet(1, -1, -1, 0, 0,0,0), 0)
    ),
    // TOwner_(lookupOrDefault<Switch>("updateT", TOwner_)),    // TOwner would be true fs the function lookupOrConstruct is called - not implemented here
    // If TOwner is zero, T will never be updated from the enthalpy or energy field in heRhoThermo or hePsiThermo
    // Since I calculate the Teqn I don't want to update T due to a calculated enthalpy
    TOwner_(lookupOrDefault<Switch>("updateT", false)),
    dpdt_(lookupOrDefault<Switch>("dpdt", true))
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// TODO: Not sure about that - Why defining it here again in this way?
Foam::autoPtr<Foam::basicThermo> Foam::basicThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return New<basicThermo>(mesh, phaseName);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicThermo::~basicThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::basicThermo& Foam::basicThermo::lookupThermo
(
    const fvPatchScalarField& pf
)
{
    // const basicThermo* thermo = pf.db().lookupObject<basicThermo>(dictName);

    // if (thermo)
    // {
    //     return *thermo;
    // }

    // Go through all thermo and if the he field and pf are on the same mesh take this thermo model
    HashTable<const basicThermo*> thermos =
        pf.db().lookupClass<basicThermo>();

    forAllConstIter(HashTable<const basicThermo*>, thermos, iter)
    {
        if
        (
            &(iter()->he().internalField())
         == &(pf.internalField())
        )
        {
            return *iter();
        }
    }

    return pf.db().lookupObject<basicThermo>(dictName);
}



void Foam::basicThermo::validate
(
    const string& app,
    const word& a
) const
{
    if (!(he().name() == phasePropertyName(a)))
    {
        FatalErrorInFunction
            << "Supported energy type is " << phasePropertyName(a)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}

void Foam::basicThermo::validate
(
    const string& app,
    const word& a,
    const word& b
) const
{
    if
    (
       !(
            he().name() == phasePropertyName(a)
         || he().name() == phasePropertyName(b)
        )
    )
    {
        FatalErrorInFunction
            << "Supported energy types are " << phasePropertyName(a)
            << " and " << phasePropertyName(b)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}

void Foam::basicThermo::validate
(
    const string& app,
    const word& a,
    const word& b,
    const word& c
) const
{
    if
    (
       !(
            he().name() == phasePropertyName(a)
         || he().name() == phasePropertyName(b)
         || he().name() == phasePropertyName(c)
        )
    )
    {
        FatalErrorInFunction
            << "Supported energy types are " << phasePropertyName(a)
            << ", " << phasePropertyName(b)
            << " and " << phasePropertyName(c)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}

void Foam::basicThermo::validate
(
    const string& app,
    const word& a,
    const word& b,
    const word& c,
    const word& d
) const
{
    if
    (
       !(
            he().name() == phasePropertyName(a)
         || he().name() == phasePropertyName(b)
         || he().name() == phasePropertyName(c)
         || he().name() == phasePropertyName(d)
        )
    )
    {
        FatalErrorInFunction
            << "Supported energy types are " << phasePropertyName(a)
            << ", " << phasePropertyName(b)
            << ", " << phasePropertyName(c)
            << " and " << phasePropertyName(d)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}

Foam::wordList Foam::basicThermo::splitThermoName
(
    const word& thermoName,
    const int nCmpt
)
{
    wordList cmpts(nCmpt);

    string::size_type beg=0, end=0, endb=0, endc=0;
    int i = 0;

    while
    (
        (endb = thermoName.find('<', beg)) != string::npos
     || (endc = thermoName.find(',', beg)) != string::npos
    )
    {
        if (endb == string::npos)
        {
            end = endc;
        }
        else if ((endc = thermoName.find(',', beg)) != string::npos)
        {
            end = std::min(endb, endc);
        }
        else
        {
            end = endb;
        }

        if (beg < end)
        {
            cmpts[i] = thermoName.substr(beg, end-beg);
            cmpts[i++].replaceAll(">","");

            // If the number of number of components in the name
            // is greater than nCmpt return an empty list
            if (i == nCmpt)
            {
                return wordList::null();
            }
        }
        beg = end + 1;
    }

    // If the number of number of components in the name is not equal to nCmpt
    // return an empty list
    if (i + 1 != nCmpt)
    {
        return wordList::null();
    }

    if (beg < thermoName.size())
    {
        cmpts[i] = thermoName.substr(beg, string::npos);
        cmpts[i].replaceAll(">","");
    }

    return cmpts;
}



Foam::volScalarField& Foam::basicThermo::p()
{
    return p_;
}


const Foam::volScalarField& Foam::basicThermo::p() const
{
    return p_;
}


const Foam::volScalarField& Foam::basicThermo::T() const
{
    return T_;
}


Foam::volScalarField& Foam::basicThermo::T()
{
    return T_;
}


const Foam::volScalarField& Foam::basicThermo::alpha() const
{
    return alpha_;
}


const Foam::scalarField& Foam::basicThermo::alpha(const label patchi) const
{
    return alpha_.boundaryField()[patchi];
}

bool Foam::basicThermo::read()
{
    return regIOobject::read();
}


// ************************************************************************* //
