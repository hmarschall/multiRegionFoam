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
#include "IOobjectList.H"


// * * * * * * * * * * * * * * * Private functions * * * * * * * * * * * * * //

template<class T>
void Foam::multiRegionSystem::assembleAndSolveBlockMatrix
(
    PtrList<GeometricField<T, fvPatchField, volMesh> >& flds,
    word fldName
) const
{
    // get number of coupled fields per field name
    label nEqns = 0;

    forAll (flds, fldI)
    {
        if(flds[fldI].name() == fldName)
        {
            nEqns++;
        }
    }

    Info << "Number of coupled fields : " << nEqns << endl;

    // assemble and solve block matrix system
    coupledFvScalarMatrix coupledEqns(nEqns);

    forAll (flds, fldI)
    {
        const fvMatrix<T>& eqn = 
            flds[fldI].mesh()
            .objectRegistry::lookupObject<IOReferencer<fvMatrix<T> > >
            (
                fldName + "Eqn"
            )();

        // Add region equation
        Info<< "Adding matrix " 
            << fldName + "Eqn" << " to blockMatrix" 
            << endl;

        // coupledEqns.set(fldI, eqn);

        coupledEqns.set(fldI, &const_cast<fvMatrix<T>& >(eqn));
    }

    coupledEqns.solve(mesh_.solutionDict().solver(fldName + "coupled"));
}

template<class T>
void Foam::multiRegionSystem::assembleAndSolveEqns
(
    PtrList<GeometricField<T, fvPatchField, volMesh> >& flds,
    regionTypeList& regions,
    word fldName
) const
{
    // assemble and solve all matrices
    forAll (regions, regI) // go through all regions
    {
        fvMatrix<T>& eqn =
            regions[regI].getCoupledEqn<T>
            (
                fldName + regions[regI].name() + "Eqn"
            );

        Info<< "Solving for " << eqn.psi().name() 
            << " in " << regions()[regI].name()
            << endl;

        eqn.relax();
        eqn.solve();
    }

//    forAll (flds, fldI)
//    {
//        if
//        (
//            flds[fldI].mesh()
//            .objectRegistry::foundObject<IOReferencer<fvMatrix<T> > >
//            (
//                fldName + "Eqn"
//            )
//        )
//        {
//        const fvMatrix<T>& eqn = 
//            flds[fldI].mesh()
//            .objectRegistry::lookupObject<IOReferencer<fvMatrix<T> > >
//            (
//                fldName + "Eqn"
//            )();

//        fvMatrix<T>& eqnRef = const_cast<fvMatrix<T>& >(eqn);

//        Info << "Solving for " << eqnRef.psi().name() << endl;

////        eqnRef.relax();
////        eqnRef.solve();
//        }
//    }
}

template<class T>
void Foam::multiRegionSystem::assembleCoupledFields
(
    PtrList<GeometricField<T, fvPatchField, volMesh> >& flds,
    regionTypeList& regions,
    const hashedWordList& fldNames
) const
{
    label nFlds = 0;

    forAll (regions, regI) // go through all regions
    {
        // get list of all objects registered to region
        IOobjectList objects
        (
            regions[regI],
            "0"
            // regions()[regI].time().timeName()
        );

        // get list of field objects of requested type
        IOobjectList volTypeObjects = 
            objects.lookupClass
            (
                GeometricField<T, fvPatchField, volMesh>::typeName
            );

        for
        (
            IOobjectList::iterator iter = volTypeObjects.begin();
            iter != volTypeObjects.end();
            ++iter
        )
        {
            // Info << "Name of coupled field : " << iter()->name() << endl;

            if
            (
                // make sure to pick up only coupled flds of requested type
                fldNames.contains(iter()->name())
            )
            {
                const GeometricField<T, fvPatchField, volMesh>& fld =
                    regions[regI].thisDb().lookupObject
                    <GeometricField<T, fvPatchField, volMesh> >
                    (
                        iter()->name()
                    );

                flds.setSize(flds.size() + 1); //no append

                flds.set(nFlds, fld);

                nFlds++;

                Info<< "Coupled field type : " 
                    << GeometricField<T, fvPatchField, volMesh>::typeName
                    << ", name : " << iter()->name()
                    << endl;
            }
        }
    }
}

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

    scalarFlds_(0),
    vectorFlds_(0),
    symmTensorFlds_(0),
    tensorFlds_(0),

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
    // Solve for regions' physics
    regions_->solveRegion();

    // Set coupled equations
    regions_->setCoupledEqns();

    // Solve region-region coupling (partitioned)

//    forAll (regions(), regI) // go through all regions
//    {
//        fvMatrix<scalar>& TEqn =
//            regions()[regI].getCoupledEqn<scalar>
//            (
//                "T" + regions()[regI].name() + "Eqn"
//            );

//        Info<< "Solving for " << TEqn.psi().name() 
//            << " in " << regions()[regI].name()
//            << endl;

//        TEqn.relax();
//        TEqn.solve();
//    }

    //- get unique list of coupled field names
    hashedWordList pcfldNames;

    forAllConstIter(fieldsTable, interfaces_->partitionedCoupledFields(), iter)
    {
        forAll(iter(), fldNameI)
        {
            word fldName = iter()[fldNameI];

            if (!pcfldNames.contains(fldName))
            {
                pcfldNames.append(fldName);
            }
        }
    }

    Info<< "List of partitioned-coupled field names : " 
        << pcfldNames << endl;

    //- assemble coupled fields by type
    scalarFlds_.clear();
    vectorFlds_.clear();
    symmTensorFlds_.clear();
    tensorFlds_.clear();

    assembleCoupledFields<scalar>(scalarFlds_, regions(), pcfldNames);
    assembleCoupledFields<vector>(vectorFlds_, regions(), pcfldNames);
    assembleCoupledFields<tensor>(tensorFlds_, regions(), pcfldNames);
    assembleCoupledFields<symmTensor>(symmTensorFlds_, regions(), pcfldNames);

    //- assemble and solve set of coupled equation (partitioned)
    forAll (pcfldNames, fldI)
    {
        word fldName = pcfldNames[fldI];

        assembleAndSolveEqns<scalar>(scalarFlds_, regions(), fldName);
//        assembleAndSolveEqns<vector>(vectorFlds_, regions(), fldName);
//        assembleAndSolveEqns<tensor>(tensorFlds_, regions(), fldName);
//        assembleAndSolveEqns<symmTensor>(symmTensorFlds_, regions(), fldName);
    }


//    // Solve region-region coupling (monolithic)
//    interfaces_->attach();

//    //- get unique list of coupled field names
//    hashedWordList mcfldNames;

//    forAllConstIter(fieldsTable, interfaces_->monolithicCoupledFields(), iter)
//    {
//        // Info << iter() << endl;

//        forAll(iter(), fldNameI)
//        {
//            word fldName = iter()[fldNameI];

//            if (!mcfldNames.contains(fldName))
//            {
//                mcfldNames.append(fldName);
//            }
//        }
//    }

//    Info<< "List of monolithic-coupled field names : " 
//        << mcfldNames << endl;

//    //- assemble coupled fields by type
//    scalarFlds_.clear();
//    vectorFlds_.clear();
//    symmTensorFlds_.clear();
//    tensorFlds_.clear();

//    assembleCoupledFields<scalar>(scalarFlds_, regions(), mcfldNames);
//    assembleCoupledFields<vector>(vectorFlds_, regions(), mcfldNames);
//    assembleCoupledFields<tensor>(tensorFlds_, regions(), mcfldNames);
//    assembleCoupledFields<symmTensor>(symmTensorFlds_, regions(), mcfldNames);

//    //- assemble and solve coupled block-matrix systems
//    forAll (mcfldNames, fldI)
//    {
//        word fldName = mcfldNames[fldI];

//        assembleAndSolveBlockMatrix<scalar>(scalarFlds_, fldName);
//        assembleAndSolveBlockMatrix<vector>(vectorFlds_, fldName);
//        assembleAndSolveBlockMatrix<tensor>(tensorFlds_, fldName);
//        assembleAndSolveBlockMatrix<symmTensor>(symmTensorFlds_, fldName);
//    }

//    // interfaces_->detach();
}

void Foam::multiRegionSystem::setCoupledEqns()
{
    regions_->setCoupledEqns();

    //- solve partitioned inter-region coupling
//    for (int coupleIter=1; coupleIter<=maxCoupleIter_; coupleIter++)
//    {

//    }
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
