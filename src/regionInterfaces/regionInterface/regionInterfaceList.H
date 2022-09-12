/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
     \\/     M anipulation  |
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

Class
    Foam::regionInterfaceList

Description
    List container for region interfaces

SourceFiles
    regionInterfaceList.C

\*---------------------------------------------------------------------------*/

#ifndef regionInterfaceList_H
#define regionInterfaceList_H

#include "Tuple2.H"
#include "fvMesh.H"
#include "dictionary.H"
#include "hashedWordList.H"
#include "regionInterfaceProperties.H"
#include "regionInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class regionInterfaceList Declaration
\*---------------------------------------------------------------------------*/

class regionInterfaceList
:
    public PtrList<regionInterface>
{
public:

    typedef
        HashTable
        <
            wordList, interfaceKey, interfaceKey::hash
        >
        fieldsTable;

private:

    label index_;

    wordList interfaceNames_;

    fieldsTable monolithicCoupledFields_;

    fieldsTable partitionedCoupledFields_;

    typedef List<Tuple2<Pair<Pair<word> >, wordList> > interfaceList;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionInterfaceList(const regionInterfaceList&);

        //- Disallow default bitwise assignment
        void operator=(const regionInterfaceList&);


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
    regionInterfaceList(const Time& runTime);

    // Destructor
    virtual ~regionInterfaceList()
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
        
        //- Interfacial transport
        void makeUs() const;
        void makeK() const;
        void makePhis() const;
        
        //- Interface properties
        void makeFaMesh() const;

        //- Get min time step
        scalar getMinDeltaT();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
