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

#include "janafThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class equationOfState>
void Foam::janafThermo<equationOfState>::checkInputData() const
{
    if (Tlow_ >= Thigh_)
    {
        FatalErrorInFunction
            << "Tlow(" << Tlow_ << ") >= Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }

    if (Tcommon_ <= Tlow_)
    {
        FatalErrorInFunction
            << "Tcommon(" << Tcommon_ << ") <= Tlow(" << Tlow_ << ')'
            << exit(FatalError);
    }

    if (Tcommon_ > Thigh_)
    {
        FatalErrorInFunction
            << "Tcommon(" << Tcommon_ << ") > Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class equationOfState>
Foam::janafThermo<equationOfState>::janafThermo(const dictionary& dict)
:
    equationOfState(dict),
    Tlow_(readScalar(dict.subDict("thermodynamics").lookup("Tlow"))),
    Thigh_(readScalar(dict.subDict("thermodynamics").lookup("Thigh"))),
    Tcommon_(readScalar(dict.subDict("thermodynamics").lookup("Tcommon"))),
    highCpCoeffs_(dict.subDict("thermodynamics").lookup("highCpCoeffs")),
    lowCpCoeffs_(dict.subDict("thermodynamics").lookup("lowCpCoeffs"))
{
    // Convert coefficients to mass-basis
    for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        highCpCoeffs_[coefLabel] *= this->R();
        lowCpCoeffs_[coefLabel] *= this->R();
    }

    checkInputData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class equationOfState>
void Foam::janafThermo<equationOfState>::write(Ostream& os) const
{
    equationOfState::write(os);

    // Convert coefficients back to dimensionless form
    coeffArray highCpCoeffs;
    coeffArray lowCpCoeffs;
    for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        highCpCoeffs[coefLabel] = highCpCoeffs_[coefLabel]/this->R();
        lowCpCoeffs[coefLabel] = lowCpCoeffs_[coefLabel]/this->R();
    }

    // // Entries in dictionary format
    // {
    //     os.beginBlock("thermodynamics");
    //     os.writeEntry("Tlow", Tlow_);
    //     os.writeEntry("Thigh", Thigh_);
    //     os.writeEntry("Tcommon", Tcommon_);
    //     os.writeEntry("highCpCoeffs", highCpCoeffs);
    //     os.writeEntry("lowCpCoeffs", lowCpCoeffs);
    //     os.endBlock();
    // }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class equationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const janafThermo<equationOfState>& jt
)
{
    os  << static_cast<const equationOfState&>(jt) << nl
        << "    " << jt.Tlow_
        << tab << jt.Thigh_
        << tab << jt.Tcommon_;

    os << nl << "    ";

    for
    (
        label coefLabel = 0;
        coefLabel<janafThermo<equationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        os << jt.highCpCoeffs_[coefLabel] << ' ';
    }

    os << nl << "    ";

    for
    (
        label coefLabel = 0;
        coefLabel<janafThermo<equationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        os << jt.lowCpCoeffs_[coefLabel] << ' ';
    }

    os << endl;

    os.check
    (
        "operator<<(Ostream& os, const janafThermo<equationOfState>& jt)"
    );

    return os;
}


// ************************************************************************* //
