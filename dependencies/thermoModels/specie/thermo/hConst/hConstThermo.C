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

#include "hConstThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class equationOfState>
Foam::hConstThermo<equationOfState>::hConstThermo(const dictionary& dict)
:
    equationOfState(dict),
    Cp_(readScalar(dict.subDict("thermodynamics").lookup("Cp"))),
    Hf_(readScalar(dict.subDict("thermodynamics").lookup("Hf")))
{
    Info << "Print Cp_ in hConstThermo.C " << Cp_ << endl;
    Info << "Print Hf_ in hConstThermo.C " << Hf_ << endl;

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class equationOfState>
void Foam::hConstThermo<equationOfState>::write(Ostream& os) const
{
    equationOfState::write(os);

    // // Entries in dictionary format
    // {
    //     os.beginBlock("thermodynamics");
    //     os.writeEntry("Cp", Cp_);
    //     os.writeEntry("Hf", Hf_);
    //     os.endBlock();
    // }
}

// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class equationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const hConstThermo<equationOfState>& ct
)
{
    os  << static_cast<const equationOfState&>(ct) << tab
        << ct.Cp_ << tab << ct.Hf_;

    os.check("Ostream& operator<<(Ostream& os, const hConstThermo& ct)");
    return os;
}


// ************************************************************************* //
