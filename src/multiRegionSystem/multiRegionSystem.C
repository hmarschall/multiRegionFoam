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

#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::multiRegionSystem::couplingMethods,
        2
    >::names[] =
    {
        "partitioned",
        "coupled"
    };
}

const Foam::NamedEnum
<
    Foam::multiRegionSystem::couplingMethods,
    2
> Foam::multiRegionSystem::couplingMethodsNames_;


// * * * * * * * * * * * * * * * Private functions * * * * * * * * * * * * * //

template< template<class> class M, class T>
void Foam::multiRegionSystem::assembleAndSolveCoupledMatrix
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

    if (nEqns == 0) return;

    Info<< fldName << ": "
       << "Number of coupled fields : " << nEqns << endl;

    // assemble and solve block matrix system
    coupledFvMatrix<T> coupledEqns(nEqns);

    // assemble and solve all matrices one-by-one
    label nReg = 0;
    forAll (regions_(), regI)
    {
        // Check if coupled field is registered to region 
        // and if it is of correct type
        if
        (
            !(
                regions_()[regI].mesh().thisDb().foundObject
                <GeometricField<T, fvPatchField, volMesh> >
                (
                    fldName
                )
            )
        )
        {
            continue;
        }

        regionType& rg = const_cast<regionType&>(regions_()[regI]);

        M<T>* eqn =
            rg.getCoupledEqn<M,T>
            (
                fldName + rg.mesh().name() + "Eqn"
            );

        coupledEqns.set(nReg, eqn);

        nReg++;
    }

    coupledEqns.solve
    (
        regions_()[0].mesh().solutionDict().solver(fldName + "coupled")
    );
}

template< template<class> class M, class T>
void Foam::multiRegionSystem::assembleAndSolveEqns
(
    word fldName
) const
{
    // assemble and solve all matrices one-by-one
    forAll (regions_(), regI)
    {
        regionType& rg = const_cast<regionType&>(regions_()[regI]);
        // Check if coupled field is registered to region mesh
        // and if it is of correct type
        // and if this region holds the equation
        if
        (
            !(
                rg.foundCoupledEqn
                (
                    fldName + rg.mesh().name() + "Eqn"
                ) 
                && 
                rg.mesh().thisDb().foundObject
                <GeometricField<T, fvPatchField, volMesh> >
                (
                    fldName
                )
            )
        )
        {
            continue;
        }

        // It is an error to use regionCouplePolyPatch
        // TODO: need more consistency checks:
        //  - Is the regionCouple bc set for fldName?
        //  - AND: Is this region adjacent to the relevant regionInterface?
        dynamicFvMesh& mesh = const_cast<dynamicFvMesh&>(regions_()[regI].mesh());

        {
            const polyPatchList& patches = mesh.boundaryMesh();

            forAll (patches, patchI)
            {
                if (isType<regionCouplePolyPatch>(patches[patchI]))
                {
                    FatalError  << "Error: Attempt to solve partitioned " 
                                << "with a patch of type " 
                                << regionCouplePolyPatch::typeName
                                << exit(FatalError);
                }
            }

            // Force recalculation of weights
            mesh.surfaceInterpolation::movePoints();
        }

        // set and get coupled equation with the possibly updated boundary conditions
        rg.setCoupledEqns();
        auto* eqn =
            rg.getCoupledEqn<M,T>
            (
                fldName + rg.mesh().name() + "Eqn"
            );

        Info<< nl 
            << "Solving for " << eqn->psi().name() 
            << " in " << rg.mesh().name()
            << endl;


        pimpleControl pimpleControlRegion(rg.mesh());

        while (pimpleControlRegion.correctNonOrthogonal())
        {
            rg.relaxEqn<T>(eqn);

            eqn->solve();

            //rg.postSolve();

            // Assemble and get equation for new subloop step
            rg.setCoupledEqns();
            eqn = rg.getCoupledEqn<M,T>(fldName + rg.mesh().name() + "Eqn");
        }
    }
}

template<class T>
void Foam::multiRegionSystem::assembleCoupledFields
(
    List<PtrList<GeometricField<T, fvPatchField, volMesh> > >& flds,
    const List<hashedWordList>& fldNms
) const
{
//    label npcFlds = 0;
//    label nmcFlds = 0;
    label n = 0;

    forAll (couplingMethodsNames_, cplI)
    {
        n = 0;

        forAll (regions_(), regI) // go through all regions
        {
            // get list of all objects registered to region
            IOobjectList objects
            (
                regions_()[regI].mesh(),
                "0"
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
                if
                (
                    fldNms[cplI].contains(iter()->name())
                )
                {
                    const GeometricField<T, fvPatchField, volMesh>& fld =
                        regions_()[regI].mesh().thisDb().lookupObject
                        <GeometricField<T, fvPatchField, volMesh> >
                        (
                            iter()->name()
                        );

                    flds[cplI].setSize(flds[cplI].size() + 1);

                    flds[cplI].set(n, fld);

                    n++;

                    // Info<< "Name of coupled field in region "
                    //    << regions_()[regI].name() << " : " 
                    //    << iter()->name() << endl;
                }
            }
        }
    } // end cpl methods
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiRegionSystem::multiRegionSystem
(
    const Time& runTime
)
:
    IOdictionary
    (
        IOobject
        (
            "multiRegionProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    runTime_(runTime),

    regions_(),
    interfaces_(),

    scalarFlds_(2),
    vectorFlds_(2),
    symmTensorFlds_(2),
    tensorFlds_(2),
    vector4Flds_(2),

    fldNames_(2),

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
            runTime_
        )
    );

    Info << nl << "Creating region interfaces:" << endl;
    interfaces_.set
    (
        new regionInterfaceList
        (
            runTime_
        )
    );

    //- set unique list of coupled field names
    //  (0 = partitioned, 1 = monolithic)
    fldNames_[0] = interfaces_->pcFldNames();
    fldNames_[1] = interfaces_->mcFldNames();

    //- assemble list of coupled fields
    assembleCoupledFields<scalar>(scalarFlds_, fldNames_);
    assembleCoupledFields<vector>(vectorFlds_, fldNames_);
    assembleCoupledFields<tensor>(tensorFlds_, fldNames_);
    assembleCoupledFields<symmTensor>(symmTensorFlds_, fldNames_);
    assembleCoupledFields<vector4>(vector4Flds_, fldNames_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiRegionSystem::~multiRegionSystem()
{
    regions_->clear();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::multiRegionSystem::preSolve()
{
    // Correct region properties and update meshes
    regions_->preSolve();

    // Update interfaces on mesh change (motion or topology)
    interfaces_->update();
}


void Foam::multiRegionSystem::setRDeltaT()
{
    regions_->setRDeltaT();
}


void Foam::multiRegionSystem::solve()
{
    // Set coupled equations
    regions_->setCoupledEqns();

    interfaces_->detach();

    // Solve each region physics
    regions_->solveRegion();

    // Solve pressure-velocity system using PIMPLE
    regions_->solvePIMPLE();

    // Solve region-region coupling (partitioned)
    forAll (fldNames_[0], fldI)
    {
        word fldName = fldNames_[0][fldI];

        // outer coupling loop
        for (int coupleIter=1; coupleIter<=maxCoupleIter_; coupleIter++)
        {
            Info << "(1)" << endl;
            assembleAndSolveEqns<fvMatrix, scalar>(fldName);
            Info << "(2)" << endl;
            assembleAndSolveEqns<fvMatrix, vector>(fldName);
            assembleAndSolveEqns<fvMatrix, tensor>(fldName);
            assembleAndSolveEqns<fvBlockMatrix, vector4>(fldName);
            // assembleAndSolveEqns<symmTensor>(fldName);
        }
    }


    // Solve region-region coupling (monolithic)
    interfaces_->attach();

    forAll (fldNames_[1], fldI)
    {
        word fldName = fldNames_[1][fldI];

        assembleAndSolveCoupledMatrix<fvMatrix, scalar>(scalarFlds_[1], fldName);
        assembleAndSolveCoupledMatrix<fvMatrix, vector>(vectorFlds_[1], fldName);
        assembleAndSolveCoupledMatrix<fvMatrix, tensor>(tensorFlds_[1], fldName);
        // assembleAndSolveCoupledMatrix<symmTensor>(symmTensorFlds_[1], fldName);
        // assembleAndSolveCoupledMatrix<fvBlockMatrix, vector4>(vector4Flds_[1], fldName);
    }
}

void Foam::multiRegionSystem::setCoupledEqns()
{
    regions_->setCoupledEqns();
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
