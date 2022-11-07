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
    monolithicCoupledFields_(),
    partitionedCoupledFields_(),
    runTime_(runTime),
    monolithicTypeInterfaces_(runTime_.time(), "monolithic"),
    partitionedTypeInterfaces_(runTime_.time(), "partitioned"),
    pcFldNames_(),
    mcFldNames_()
{
    this->setSize
    (
        partitionedTypeInterfaces_.size()
      + monolithicTypeInterfaces_.size()
    );

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
    forAll (rip, cpldPatchI)
    {
        const dictionary& dict = rip[cpldPatchI].dict();

        word interfaceType(dict.lookup("interfaceType"));
        coupledFields fields(dict.lookup("coupledFields"));
        coupledPatchPair patchPair(dict.lookup("coupledPatchPair"));
        dictionary interfaceDict(dict.subDict(interfaceType + "Coeffs"));

        if (patchPair.size() != 2)
        {
            FatalErrorIn
            (
                "regionInterfaceList::reset"
            )   << "An interface is made of two patches." << nl
                << "For each interface create a separate entry."
                << abort(FatalError);
        }

        // first patch
        const fvMesh& firstRegion = 
            runTime_.lookupObject<fvMesh>
            (
                patchPair[0].first()
            );

        label firstPatchID = 
            firstRegion.boundaryMesh().findPatchID
            (
                patchPair[0].second()
            );

        const fvPatch& firstPatch = 
            firstRegion.boundary()[firstPatchID];

        // second patch
        const fvMesh& secondRegion = 
            runTime_.lookupObject<fvMesh>
            (
                patchPair[1].first()
            );

        label secondPatchID = 
            secondRegion.boundaryMesh().findPatchID
            (
                patchPair[1].second()
            );

        const fvPatch& secondPatch = 
            secondRegion.boundary()[secondPatchID];

        this->set
        (
            index_++,
            regionInterface::New
            (
                interfaceType,
                interfaceDict,
                runTime_,
                firstPatch,
                secondPatch
            )
        );
    }
}


void Foam::regionInterfaceList::setFieldNamesPartitionedCoupling
(
    const regionInterfaceProperties& rip
)
{
    forAll (rip, cpldPatchI)
    {
        const dictionary& dict = rip[cpldPatchI].dict();

        coupledPatchPair patchPair(dict.lookup("coupledPatchPair"));
        coupledFields fields(dict.lookup("coupledFields"));

        const interfaceKey key
        (
            patchPair[0].first() + patchPair[0].second(),
            patchPair[1].first() + patchPair[1].second()
        );

        partitionedCoupledFields_.insert
        (
            key,
            fields
        );
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
    forAll (rip, cpldPatchI)
    {
        const dictionary& dict = rip[cpldPatchI].dict();

        coupledPatchPair patchPair(dict.lookup("coupledPatchPair"));
        coupledFields fields(dict.lookup("coupledFields"));

        const interfaceKey key
        (
            patchPair[0].first() + patchPair[0].second(),
            patchPair[1].first() + patchPair[1].second()
        );

        monolithicCoupledFields_.insert
        (
            key,
            fields
        );
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
   forAll(*this, i)
   {
       if(this->operator[](i).changing())
       {
           this->operator[](i).updateInterpolatorAndGlobalPatches();
       }
   }
}

Foam::scalar Foam::regionInterfaceList::getMinDeltaT()
{
    scalar minDeltaT = GREAT;
    forAll(*this, i)
    {
        minDeltaT = min(minDeltaT, this->operator[](i).getMinDeltaT());
    }

    return minDeltaT;
}


// ************************************************************************* //
