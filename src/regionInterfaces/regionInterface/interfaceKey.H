/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Fork:     foam-extend
    \\  /    A nd           | Version:  4.1                                 
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
    Foam::interfaceKey

SourceFiles

Author
    Holger Marschall (holger.marschall@tu-darmstadt.de, Affiliation A)
    Constantin Habes (constantin.habes@tu-darmstadt.de, Affiliation A)
    Heba Alkafri (heba.alkafri@tu-darmstadt.de, Affiliation A)

Contact
    Holger Marschall (holger.marschall@tu-darmstadt.de)
    main developer and principal investigator
    TU Darmstadt

Affiliations
    Affiliation A)
    Computational Multiphase Flow
    Department of Mathematics
    Technical University of Darmstadt, Germany

Acknowledgement
    Funded by
    Hessian Ministry of Higher Education, Research, Science and the Arts,
    National High Performance Computing Center for Computational Engineering
    Science (NHR4CES)

Description

    This file is part of the multiRegionFoam library.

    You may refer to this software as :
    https://doi.org/10.48550/arXiv.2306.01924

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#ifndef interfaceKey_H
#define interfaceKey_H

#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class interfaceKey;

bool operator==(const interfaceKey& a, const interfaceKey& b);
bool operator!=(const interfaceKey& a, const interfaceKey& b);

Istream& operator>>(Istream& is, interfaceKey& key);
Ostream& operator<<(Ostream& os, const interfaceKey& key);

/*---------------------------------------------------------------------------*\
                         Class interfaceKey Declaration
\*---------------------------------------------------------------------------*/

class interfaceKey
:
    public Pair<word>
{
    // Private data

        //- Flag to indicate whether ordering is important
        bool ordered_;

    public:

        //- Ordered or unordered hashing of word pair
        struct hash
        {
            //- Generate a hash from a phase pair key
            label operator()(const interfaceKey& key) const;
        };

    // Constructors

        //- Construct null
        interfaceKey() {} // = default

        //- Construct from names and optional ordering flag
        interfaceKey
        (
            const word& name1,
            const word& name2,
            const bool ordered = false
        );

    //- Destructor
//    virtual ~interfaceKey() = default;
    virtual ~interfaceKey()
    {}

    // Access

        //- Return the ordered flag
        bool ordered() const;

    // Friend Operators

        //- Test for equality
        friend bool operator==(const interfaceKey& a, const interfaceKey& b);

        //- Test for inequality
        friend bool operator!=(const interfaceKey& a, const interfaceKey& b);

        //- Read from Istream
        friend Istream& operator>>(Istream& is, interfaceKey& key);

        //- Write to Ostream
        friend Ostream& operator<<(Ostream& os, const interfaceKey& key);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
