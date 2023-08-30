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

Description
    Base class of the thermophysical property types.

\*---------------------------------------------------------------------------*/

#include "specie.H"
#include "IOstreams.H"
#include "dimensionedConstants.H"

/* * * * * * * * * * * * * public constants  * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(specie, 0);
}

//- Universal gas constant (default in [J/(kmol K)])
//- kmol is applied here to compensate for the calculation of
// Ri = R/M , where M is given in g/mol
const Foam::scalar
Foam::specie::RR
(
    8314.51
    // "Universal gas constant [J/(kmol K)]"
);


//- Standard pressure (default in [Pa])
const Foam::scalar
Foam::specie::Pstd
(
    // "Pstd",
    1.0e5
    // "Standard pressure [Pa]"
);

//- Standard temperature (default in [K])
const Foam::scalar
Foam::specie::Tstd
(
    // "Tstd",
    298.15
    // "Standard temperature [K]"
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::specie::specie(const dictionary& dict)
:
    name_(dict.dictName()),
    Y_(dict.subDict("specie").lookupOrDefault<scalar>("Y", 1)),
    molWeight_(readScalar(dict.subDict("specie").lookup("molWeight"))),
    molarVolume_(readScalar(dict.subDict("specie").lookup("molarVolume"))),
    z_(readScalar(dict.subDict("specie").lookup("valenceNumber")))
{
    // NOTE:
    // E.g. the name here in FE is with this 
    Info << "Print constructor specie" << endl;
    Info << "Print constructor specie: name_ " << name_ << endl;
    // Print constructor specie: name_ thermophysicalProperties::O2      
    // Kann in vsCode als String auch nicht angezeigt werden...
    Info << "Debug line " << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::specie::write(Ostream& os) const
{
    // // Entries in dictionary format
    // {
    //     os.beginBlock("specie");
    //     os.writeEntryIfDifferent<scalar>("massFraction", 1, Y_);
    //     os.writeEntry("molWeight", molWeight_);
    //     os.endBlock();
    // }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const specie& st)
{
    os  << st.name_ << tab
        << st.Y_ << tab
        << st.molWeight_;

    os.check("Ostream& operator<<(Ostream& os, const specie& st)");
    return os;
}


// ************************************************************************* //
