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

#include "hPolynomialThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class equationOfState, int PolySize>
Foam::hPolynomialThermo<equationOfState, PolySize>::hPolynomialThermo
(
    const dictionary& dict
)
:
    equationOfState(dict),
    Hf_(readScalar(dict.subDict("thermodynamics").lookup("Hf"))),
    Sf_(readScalar(dict.subDict("thermodynamics").lookup("Sf"))),
    CpCoeffs_(dict.subDict("thermodynamics").lookup(coeffsName("Cp"))),
    hCoeffs_(),
    sCoeffs_()
{
    hCoeffs_ = CpCoeffs_.integral();
    sCoeffs_ = CpCoeffs_.integralMinus1();

    // Offset h poly so that it is relative to the enthalpy at Tstd
    hCoeffs_[0] += Hf_ - hCoeffs_.value(this->Tstd);

    // Offset s poly so that it is relative to the entropy at Tstd
    sCoeffs_[0] += Sf_ - sCoeffs_.value(this->Tstd);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class equationOfState, int PolySize>
void Foam::hPolynomialThermo<equationOfState, PolySize>::write
(
    Ostream& os
) const
{
    equationOfState::write(os);

    // // Entries in dictionary format
    // {
    //     os.beginBlock("thermodynamics");
    //     os.writeEntry("Hf", Hf_);
    //     os.writeEntry("Sf", Sf_);
    //     os.writeEntry(coeffsName("Cp"), CpCoeffs_);
    //     os.endBlock();
    // }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class equationOfState, int PolySize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const hPolynomialThermo<equationOfState, PolySize>& pt
)
{
    os  << static_cast<const equationOfState&>(pt) << tab
        << pt.Hf_/pt.W() << tab
        << pt.Sf_ << tab
        << "cpPolynomial" << tab << pt.cpPolynomial_/pt.W();

    os.check
    (
        "operator<<"
        "("
            "Ostream&, "
            "const hPolynomialThermo<equationOfState, PolySize>&"
        ")"
    );

    return os;
}


// ************************************************************************* //
