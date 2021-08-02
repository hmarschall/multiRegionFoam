/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
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

#include "fvCFD.H"
#include "coupledFvMatrices.H"

#include "multiRegionSystem.H"
#include "HashPtrTable.H"
#include "regionTypeList.H"
#include "regionInterfaceList.H"

// * * * * * * * * * * * * * * * Private functions * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiRegionSystem::multiRegionSystem
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "multiRegionProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    mesh_(mesh),

    regions_(),
    interfaces_(),

    maxCoupleIter_
    (
        this->subDict("partitionedCoupling")
        .lookupOrDefault<label>("maxCoupleIter", 50)
    )
{
    Info << "Creating regions:" << endl;
    regions_.set
    (
        new regionTypeList
        (
            mesh_
        )
    );

    Info << nl << "Creating region interfaces:" << endl;
    interfaces_.set
    (
        new regionInterfaceList
        (
            mesh_
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiRegionSystem::~multiRegionSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::multiRegionSystem::correct()
{
    //- correct regions' properties
    regions_->correct();
}


void Foam::multiRegionSystem::setRDeltaT()
{
    regions_->setRDeltaT();
}


void Foam::multiRegionSystem::solve()
{
    //- solve for regions' physics
    regions_->solveRegion();

    //- solve region-region coupled (monolithic)

    // Info << interfaces_->partitionedCoupledFields() << endl;

    typedef
        HashTable
        <
            wordList, interfaceKey, interfaceKey::hash
        >
        fieldsTable;

    interfaces_->attach();

    // solve for monolithically coupled transport equations

    //- determine size
    label NEqns = 2;
    coupledFvScalarMatrix scalarEqns(NEqns);

    // word fldName = "T";

    //- assemble matrix systems

    forAllConstIter(fieldsTable, interfaces_->partitionedCoupledFields(), iter)
    {
        forAll(iter(), cfldI)
        {
            // Info << iter() << endl;

            forAll(iter(), fldNameI)
            {
                word fldName = iter()[fldNameI];

                forAll(regions(), regI)
                {
                    if (regions()[regI].foundObject<volScalarField>(fldName))
                    {
                        Info<< "Found field " << fldName
                            << " in " << regions()[regI].name()
                            << " for monolithic coupling"
                            << endl;
                        // scalarEqns.set
                        // (
                        //     regI, 
                        //     regions()[regI].coupledScalarEquations(fldName)
                        // );
                    }
                    // TODO: Sanity check
                    // else
                    // {
                    //     forAll(interfaces(), intI)
                    //     {
                    //         FatalError << "Coupling field " << fldName
                    //             << " not found"
                    //             << exit(FatalError);
                    //     }
                    // }
                }

                // scalarEqns.solve
                // (
                //     mesh_.solutionDict().solver(fldNameI + "coupled")
                // );
            }
        }
    }

    // interfaces_->detach();
}

void Foam::multiRegionSystem::solveCoupledPartitioned()
{
    //- solve partitioned inter-region coupling
    for (int coupleIter=1; coupleIter<=maxCoupleIter_; coupleIter++)
    {
        regions_->solveCoupledPartitioned();
    }
}


Foam::regionTypeList& Foam::multiRegionSystem::regions()
{
    return regions_();
}

Foam::regionInterfaceList& Foam::multiRegionSystem::interfaces()
{
    return interfaces_();
}


//bool multiRegionSystem::read()
//{
//    if (IOdictionary::read())
//    {
//        bool readOK = true;

//        lookup("coupledRegionInterfaces") >> coupledInterfaces_;

//        return readOK;
//    }

//    return false;
//}

// ************************************************************************* //
