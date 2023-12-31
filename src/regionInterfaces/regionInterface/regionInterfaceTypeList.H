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
    Foam::regionInterfaceTypeList

SourceFiles
    regionInterfaceTypeList.C

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
    List container for region interfaces

    This file is part of the multiRegionFoam library.

    You may refer to this software as :
    https://doi.org/10.48550/arXiv.2306.01924

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#ifndef regionInterfaceTypeList_H
#define regionInterfaceTypeList_H

#include "Tuple2.H"
#include "fvMesh.H"
#include "dictionary.H"
#include "hashedWordList.H"
#include "regionInterfaceProperties.H"
#include "regionInterfaceType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class regionInterfaceTypeList Declaration
\*---------------------------------------------------------------------------*/

class regionInterfaceTypeList
:
    public PtrList<regionInterfaceType>
{
public:

    typedef
        HashTable
        <
            wordList, interfaceKey, interfaceKey::hash
        >
        fieldsTable;

private:

    // Private data

        label index_;

        fieldsTable monolithicCoupledFields_;

        fieldsTable partitionedCoupledFields_;

        typedef List<Pair<word> > coupledPatchPair;

        typedef List<word> coupledFields;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionInterfaceTypeList(const regionInterfaceTypeList&);

        //- Disallow default bitwise assignment
        void operator=(const regionInterfaceTypeList&);

protected:

    // Protected data

        //- Reference to the database
        const Time& runTime_;

        //- Tables of coupling interfaces
        regionInterfaceProperties monolithicTypeInterfaces_;

        regionInterfaceProperties partitionedTypeInterfaces_;

        //- Unique list of coupled field names
        //  (p = partitioned, m = monolithic)
        hashedWordList pcFldNames_;

        hashedWordList mcFldNames_;

public:

    // Constructor
    regionInterfaceTypeList(const Time& runTime);

    // Destructor
    virtual ~regionInterfaceTypeList()
    {}

    // Member Functions

        //- Access field tables for monolithic coupling
        const fieldsTable& monolithicCoupledFields()
        {
            return monolithicCoupledFields_;
        }

        //- Access field tables for partitioned coupling
        const fieldsTable& partitionedCoupledFields()
        {
            return partitionedCoupledFields_;
        }

        //- Access coupled field names
        const hashedWordList& pcFldNames()
        {
            return pcFldNames_;
        }

        const hashedWordList& mcFldNames()
        {
            return mcFldNames_;
        }

        //- Return active status
//        bool coupled(const bool coupled = false) const;

        //- Reset the region interface list
        void reset
        (
            const regionInterfaceProperties& regionInterfaces
        );

        //- Set field names to couple in partitioned fashion
        void setFieldNamesPartitionedCoupling
        (
            const regionInterfaceProperties& regionInterfaces
        );

        //- Set field names to couple in monolithic fashion
        void setFieldNamesMonolithicCoupling
        (
            const regionInterfaceProperties& regionInterfaces
        );

        //- Attach interface
        void attach();

        //- Detach interface
        void detach();

        //- Update boundary
        void update();

        //- Get min time step over regions
        scalar getMinDeltaT();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
