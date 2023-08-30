/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
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

#include "rhoConst.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::rhoConst<Specie>::rhoConst(const dictionary& dict)
:
    Specie(dict),
    rho_(readScalar(dict.subDict("equationOfState").lookup("rho")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::rhoConst<Specie>::write(Ostream& os) const
{
    Specie::write(os);

    // // Entries in dictionary format
    // {
    //     os.beginBlock("equationOfState");
    //     os.writeEntry("rho", rho_);
    //     os.endBlock();
    // }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<(Ostream& os, const rhoConst<Specie>& ico)
{
    os  << static_cast<const Specie&>(ico);

    os.check("Ostream& operator<<(Ostream& os, const rhoConst& st)");

    return os;
}


// ************************************************************************* //
