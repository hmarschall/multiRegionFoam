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

\*---------------------------------------------------------------------------*/

#include "regionInterfaceList.H"
#include "volFields.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionInterfaceList::regionInterfaceList
(
    const Time& runTime
)
:
    PtrList<regionInterface>(),
    index_(0),
    interfaceNames_(),
    monolithicCoupledFields_(),
    partitionedCoupledFields_(),
    runTime_(runTime),
    monolithicTypeInterfaces_(runTime_.time(), "monolithic"),
    partitionedTypeInterfaces_(runTime_.time(), "partitioned"),
    pcFldNames_(),
    mcFldNames_()
{
    if (partitionedTypeInterfaces_.size() > 0)
    {
        reset(partitionedTypeInterfaces_);

        setFieldNamesPartitionedCoupling(partitionedTypeInterfaces_);
    }

    if (monolithicTypeInterfaces_.size() > 0)
    {
        reset(monolithicTypeInterfaces_);
        setFieldNamesMonolithicCoupling(monolithicTypeInterfaces_);
    }

    // coupled(true);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//bool Foam::regionInterfaceList::coupled(const bool warn) const
//{
//    bool a = false;
//    forAll(*this, i)
//    {
//        a = a || this->operator[](i).coupled();
//    }

//    if (warn && this->size() && !a)
//    {
//        Info<< "No interfaces coupled" << endl;
//    }

//    return a;
//}

void Foam::regionInterfaceList::reset(const regionInterfaceProperties& rip)
{
    interfaceNames_.clear();

    forAllConstIter(HashTable<interfaceList>, rip, iter)
    {
        const interfaceList& interfaces = iter();

        forAll(interfaces, interfaceI)
        {
            if (findIndex(interfaceNames_, iter.key()))
            {
                interfaceNames_.setSize(interfaceNames_.size()+1);
                interfaceNames_[interfaceI] = iter.key();
            }
        }
    }

    this->setSize(interfaceNames_.size());

    forAllConstIter(HashTable<interfaceList>, rip, iter)
    {
        const interfaceList& interfaces = iter();

        Info << "Creating " << iter.key() << endl;

        forAll(interfaces, interfaceI)
        {
            Pair<word> firstRegionPatchPair = 
                interfaces[interfaceI].first().first();

            Pair<word> secondRegionPatchPair = 
                interfaces[interfaceI].first().second();

            const fvMesh& firstRegion = 
                runTime_.lookupObject<fvMesh>
                (
                    firstRegionPatchPair.first()
                );

            label firstPatchID = 
                firstRegion.boundaryMesh().findPatchID
                (
                    firstRegionPatchPair.second()
                );

            const fvPatch& firstPatch = 
                firstRegion.boundary()[firstPatchID];

            const fvMesh& secondRegion = 
                runTime_.lookupObject<fvMesh>
                (
                    secondRegionPatchPair.first()
                );

            label secondPatchID = 
                secondRegion.boundaryMesh().findPatchID
                (
                    secondRegionPatchPair.second()
                );

            const fvPatch& secondPatch = 
                secondRegion.boundary()[secondPatchID];

            this->set
            (
                index_++,
                regionInterface::New
                (
                    runTime_,
                    firstPatch,
                    secondPatch
                )
            );
        }
    }
}


void Foam::regionInterfaceList::setFieldNamesPartitionedCoupling
(
    const regionInterfaceProperties& rip
)
{
    forAllConstIter(HashTable<interfaceList>, rip, iter)
    {
        const interfaceList& interfaces = iter();

        forAll(interfaces, interfaceI)
        {
            Pair<word> firstRegionPatchPair = 
                interfaces[interfaceI].first().first();

            Pair<word> secondRegionPatchPair = 
                interfaces[interfaceI].first().second();

            word firstName2(firstRegionPatchPair.second());
            firstName2[0] = toupper(firstName2[0]);

            word secondName1(secondRegionPatchPair.first());
            secondName1[0] = toupper(secondName1[0]);

            word secondName2(secondRegionPatchPair.second());
            secondName2[0] = toupper(secondName2[0]);

            const interfaceKey key
            (
                firstRegionPatchPair.first() + firstName2,
                secondName1 + secondName2
            );

            partitionedCoupledFields_.insert
            (
                key,
                interfaces[interfaceI].second()
            );
        }
    }

    //- get unique list of coupled field names (partitioned)
    forAllConstIter(fieldsTable, partitionedCoupledFields(), iter)
    {
        forAll(iter(), fldNameI)
        {
            word fldName = iter()[fldNameI];

            if (!pcFldNames_.contains(fldName))
            {
                pcFldNames_.append(fldName);
            }
        }
    }
}

void Foam::regionInterfaceList::setFieldNamesMonolithicCoupling
(
    const regionInterfaceProperties& rip
)
{
    forAllConstIter(HashTable<interfaceList>, rip, iter)
    {
        const interfaceList& interfaces = iter();

        forAll(interfaces, interfaceI)
        {
            Pair<word> firstRegionPatchPair = 
                interfaces[interfaceI].first().first();

            Pair<word> secondRegionPatchPair = 
                interfaces[interfaceI].first().second();

            const interfaceKey key
            (
                firstRegionPatchPair.first() + firstRegionPatchPair.second(),
                secondRegionPatchPair.first() + secondRegionPatchPair.second()
            );

            monolithicCoupledFields_.insert
            (
                key,
                interfaces[interfaceI].second()
            );
        }
    }

    //- get unique list of coupled field names (monolithic)
    forAllConstIter(fieldsTable, monolithicCoupledFields(), iter)
    {
        forAll(iter(), fldNameI)
        {
            word fldName = iter()[fldNameI];

            if (!mcFldNames_.contains(fldName))
            {
                mcFldNames_.append(fldName);
            }
        }
    }

    // Info << monolithicCoupledFields_.toc() << endl;
}

void Foam::regionInterfaceList::attach()
{
   forAll(*this, i)
   {
       this->operator[](i).attach();
   }
}

void Foam::regionInterfaceList::detach()
{
   forAll(*this, i)
   {
       this->operator[](i).detach();
   }
}

void Foam::regionInterfaceList::update()
{
//    forAll(*this, i)
//    {
//        if(this->operator[](i).changing())
//        {
//            this->operator[](i).updateInterpolatorAndGlobalPatches();
//        }
//    }
}


// ************************************************************************* //
