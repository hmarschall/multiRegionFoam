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

#include "sutherlandTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class thermo>
Foam::sutherlandTransport<thermo>::sutherlandTransport(const dictionary& dict)
:
    thermo(dict),
    As(readScalar(dict.subDict("transport").lookup("As"))),
    Ts(readScalar(dict.subDict("transport").lookup("Ts")))
{}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// template<class thermo>
// Ostream& operator<<(Ostream& os, const sutherlandTransport<thermo>& st)
// {
//     os << static_cast<const thermo&>(st) << tab << st.As << tab << st.Ts;

//     os.check
//     (
//         "Ostream& operator<<(Ostream&, const sutherlandTransport<thermo>&)"
//     );

//     return os;
// }

template<class thermo>
void Foam::sutherlandTransport<thermo>::write(Ostream& os) const
{
    // os.beginBlock(this->specie::name());

    thermo::write(os);

    // // Entries in dictionary format
    // {
    //     os.beginBlock("transport");
    //     os.writeEntry("As", As_);
    //     os.writeEntry("Ts", Ts_);
    //     os.endBlock();
    // }

    // os.endBlock();
}


template<class thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const sutherlandTransport<thermo>& st
)
{
    st.write(os);
    return os;
}


// ************************************************************************* //
