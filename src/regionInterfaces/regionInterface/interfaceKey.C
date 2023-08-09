/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Fork:     foam-extend
    \\  /    A nd           | Version:  4.0
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    interface

SourceFiles
    interface.C

Author
    Holger Marschall (holger.marschall@tu-darmstadt.de, Affiliation A)

Contact
    Holger Marschall (holger.marschall@tu-darmstadt.de)
    main developer and principal investigator
    TU Darmstadt
    ---
    Martin Wörner (martin.woerner@kit.edu)
    principal investigator
    Karlsruhe Institute of Technology

Affiliations
    Affiliation A)
    Computational Multiphase Flow
    Research Group Leader
    Department of Mathematics
    Technical University of Darmstadt, Germany

Acknowledgement
    Financed by German Research Foundation (DFG) via
    SFB TRR 150, Project 237267381 - TRR 150
    SFB 1194, Procect 265191195 - SFB 1194

Description

    This file is part of the phaseFieldFoam library.

    You may refer to this software as :
    https://doi.org/10.1002/ceat.201500089
    https://doi.org/10.1016/j.cattod.2016.03.053
    https://doi.org/10.1016/j.cpc.2018.10.015

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

    This file is part of the phaseFieldFoam library.

    You may refer to this software as :
    https://doi.org/10.1002/ceat.201500089
    https://doi.org/10.1016/j.cattod.2016.03.053
    https://doi.org/10.1016/j.cpc.2018.10.015

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#include "interfaceKey.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceKey::interfaceKey
(
    const word& name1,
    const word& name2,
    const bool ordered
)
:
    Pair<word>(name1, name2),
    ordered_(ordered)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::interfaceKey::ordered() const
{
    return ordered_;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

Foam::label Foam::interfaceKey::hash::operator()
(
    const interfaceKey& key
) const
{
//    if (key.ordered_)
//    {
//        return
//            word::hash()
//            (
//                key.first(),
//                word::hash()(key.second())
//            );
//    }

    return word::hash()(key.first()) + word::hash()(key.second());
}

// * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * * //

bool Foam::operator==
(
    const interfaceKey& a,
    const interfaceKey& b
)
{
    return
    (
        ((a.first() == b.first()) && (a.second() == b.second()))
     || ((a.first() == b.second()) && (a.second() == b.first()))
    );
}

bool Foam::operator!=
(
    const interfaceKey& a,
    const interfaceKey& b
)
{
    return !(a == b);
}

// * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, interfaceKey& key)
{
    const FixedList<word, 2> temp(is);

    key.first() = temp[0];

    key.second() = temp[1];

    return is;
}

// * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const interfaceKey& key)
{
    os
//        << token::BEGIN_LIST
        << key.first()
//        << token::SPACE
//        << (key.ordered_ ? "to" : "and")
//        << token::SPACE
        << key.second();
//        << token::END_LIST;

    return os;
}

// ************************************************************************* //
