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

    //Info<< fldName << ": "
    //   << "Number of coupled fields : " << nEqns << endl;

    // assemble and solve block matrix system
    coupledFvMatrix<T> coupledEqns(nEqns);

    // assemble and solve all matrices one-by-one
    label nReg = 0;
    DynamicList<M<T>* > eqns;

    forAll (regions_(), regI)
    {
        regionType& rg = const_cast<regionType&>(regions_()[regI]);

        word matrixSystemName =
        (
            fldName
          + rg.mesh().name() + "Mesh"
          + rg.regionTypeName() + "Type"
          + "Eqn"
        );

        // Check if coupled field is registered to region mesh
        // and if it is of correct type
        if
        (
            rg.mesh().thisDb().foundObject
            <GeometricField<T, fvPatchField, volMesh> >
            (
                fldName
            )
        )
        {
            // set coupled equation
            rg.setCoupledEqns();

            // Check if this region holds the equation
            if
            (
                !rg.foundCoupledEqn
                (
                    matrixSystemName
                )
            )
            {
                continue;
            }
        }
        else
        {
            continue;
        }

        Info<< "Get equation for "
            << fldName << " in " 
            << rg.regionTypeName() << endl;

        M<T>& eqn =
            rg.getCoupledEqn<M,T>(matrixSystemName);

//        coupledEqns.set(nReg, &eqn);

        eqns.append(new M<T>(eqn));

        coupledEqns.set(nReg, eqns[nReg]);

        nReg++;
    }

    coupledEqns.solve
    (
        regions_()[0].mesh().solutionDict().solver(fldName + "coupled")
    );

    forAll (regions_(), regI)
    {
        regionType& rg = const_cast<regionType&>(regions_()[regI]);

        word matrixSystemName =
        (
            fldName
          + rg.mesh().name() + "Mesh"
          + rg.regionTypeName() + "Type"
          + "Eqn"
        );

        // Check if coupled field is registered to region 
        // and if it is of correct type
        if
        (
            !(
                rg.foundCoupledEqn
                (
                    matrixSystemName
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

        // Post-solve actions
        rg.postSolve();
    }
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

        word matrixSystemName =
        (
            fldName
          + rg.mesh().name() + "Mesh"
          + rg.regionTypeName() + "Type"
          + "Eqn"
        );

        // Check if coupled field is registered to region mesh
        // and if it is of correct type
        if
        (
            rg.mesh().thisDb().foundObject
            <GeometricField<T, fvPatchField, volMesh> >
            (
                fldName
            )
        )
        {
            // set coupled equation
            rg.setCoupledEqns();

            // Check if this region holds the equation
            if
            (
                !rg.foundCoupledEqn
                (
                    matrixSystemName
                )
            )
            {
                continue;
            }
        }
        else
        {
            continue;
        }

        // It is an error to use regionCouplePolyPatch
        // TODO: need more consistency checks:
        //  - Is the regionCouple bc set for fldName?
        //  - AND: Is this region adjacent to the relevant regionInterface?
        fvMesh& mesh = const_cast<fvMesh&>(regions_()[regI].mesh());

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

        Info<< nl 
            << "Solving for " << fldName
            << " in " << rg.regionTypeName()
            << endl;

        M<T>& eqn =
            rg.getCoupledEqn<M,T>(matrixSystemName);

        rg.relaxEqn<T>(eqn);

        eqn.solve();

        rg.postSolve();

        // Memory management:
        // clear coupled equation for transport variable
//        const GeometricField<T, fvPatchField, volMesh>& fld =
//            rg.mesh().lookupObject
//            <
//                GeometricField<T, fvPatchField, volMesh>
//            >(fldName);

//        rg.clearCoupledEqn<GeometricField<T, fvPatchField, volMesh> >
//            (fld, rg.regionTypeName());
    }
}

template<class T>
void Foam::multiRegionSystem::assembleCoupledFields
(
    PtrList<GeometricField<T, fvPatchField, volMesh> >& flds,
    const hashedWordList& fldNms
) const
{
    // Get list of meshes registered
    const fvMesh& mesh = regions_()[0].mesh();
    const objectRegistry& obr = mesh.objectRegistry::parent();

    meshList meshes = obr.lookupClass<fvMesh>();

    label n = 0;

    // go through all meshes
    forAllConstIter(meshList, meshes, iter)
    {
        const fvMesh& mesh = *iter();

        // get list of all objects registered to region
        IOobjectList objects
        (
            mesh,
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
                fldNms.contains(iter()->name())
            )
            {
                const GeometricField<T, fvPatchField, volMesh>& fld =
                    mesh.thisDb().lookupObject
                    <
                        GeometricField<T, fvPatchField, volMesh>
                    >
                    (
                        iter()->name()
                    );

                flds.setSize(flds.size() + 1);

                flds.set(n, fld);

                n++;

                //Info<< "Name of coupled field in mesh "
                //   << mesh.name() << " : " 
                //   << iter()->name() << endl;
            }
        }
    }
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

    partitionedCoupledScalarFlds_(),
    partitionedCoupledVectorFlds_(),
    partitionedCoupledSymmTensorFlds_(),
    partitionedCoupledTensorFlds_(),
    partitionedCoupledVector4Flds_(),

    monolithicCoupledScalarFlds_(),
    monolithicCoupledVectorFlds_(),
    monolithicCoupledSymmTensorFlds_(),
    monolithicCoupledTensorFlds_(),
    monolithicCoupledVector4Flds_(),

    partitionedCoupledFldNames_(),
    monolithicCoupledFldNames_(),

    dnaControls_()
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

    // set unique list of coupled field names
    partitionedCoupledFldNames_ = interfaces_->pcFldNames();
    monolithicCoupledFldNames_ = interfaces_->mcFldNames();

    // set up DNA convergence control for partitioned coupling
    forAll (partitionedCoupledFldNames_, fldI)
    {
        dnaControls_.set
        (
            partitionedCoupledFldNames_[fldI],
            new dnaControl
            (
                runTime_,
                partitionedCoupledFldNames_[fldI],
                interfaces()
            )
        );
    }

    //- assemble list of coupled fields
    assembleCoupledFields<scalar>
        (partitionedCoupledScalarFlds_, partitionedCoupledFldNames_);

    assembleCoupledFields<vector>
        (partitionedCoupledVectorFlds_, partitionedCoupledFldNames_);

    assembleCoupledFields<tensor>
        (partitionedCoupledTensorFlds_, partitionedCoupledFldNames_);

    assembleCoupledFields<symmTensor>
        (partitionedCoupledSymmTensorFlds_, partitionedCoupledFldNames_);

    assembleCoupledFields<vector4>
        (partitionedCoupledVector4Flds_, partitionedCoupledFldNames_);

    assembleCoupledFields<scalar>
        (monolithicCoupledScalarFlds_, monolithicCoupledFldNames_);

    assembleCoupledFields<vector>
        (monolithicCoupledVectorFlds_, monolithicCoupledFldNames_);

    assembleCoupledFields<tensor>
        (monolithicCoupledTensorFlds_, monolithicCoupledFldNames_);

    assembleCoupledFields<symmTensor>
        (monolithicCoupledSymmTensorFlds_, monolithicCoupledFldNames_);

    assembleCoupledFields<vector4>
        (monolithicCoupledVector4Flds_, monolithicCoupledFldNames_);
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

    // ALE mesh motion corrector
    regions_->meshMotionCorrector();

    // Update interfaces on mesh change (motion or topology)
    interfaces_->update();
}


Foam::scalar Foam::multiRegionSystem::getMinDeltaT()
{
    //- set deltaT based on volumetric stability criterion
    scalar minRegionDeltaT = regions_->getMinDeltaT();
    //- set deltaT based on interface stability criterion
    scalar minInterfaceDeltaT = interfaces_->getMinDeltaT();

    return min(minRegionDeltaT, minInterfaceDeltaT);
}


void Foam::multiRegionSystem::solve()
{
    // Set coupled equations
    // regions_->setCoupledEqns();

    interfaces_->detach();

    //- Initial cpu time call
    runTime_.cpuTimeIncrement();

    // Solve individual region physics
    regions_->solveRegion();

    Info<< "Solved region specific physics in "
        << runTime_.cpuTimeIncrement() << " s." << endl;

    // Solve pressure-velocity system using PIMPLE
    // Check if at least one region implements PIMPLE loop
    if (regions_->usesPIMPLE() && !partitionedCoupledFldNames_.contains("pUPimple"))
    {
        // PIMPLE p-U-coupling
        regions_->solvePIMPLE();

        Info<< "Solved PIMPLE without coupling in "
            << runTime_.cpuTimeIncrement() << " s." << endl;
    }

    // Solve region-region coupling (partitioned)  
    forAll (partitionedCoupledFldNames_, fldI)
    {
        word fldName = partitionedCoupledFldNames_[fldI];

        //- Solve pressure-velocity system using PIMPLE
        if (fldName == "pUPimple")
        {
            while (dnaControls_[fldName]->loop())
            {
                // PIMPLE p-U-coupling
                regions_->solvePIMPLE();

                // ALE mesh motion corrector
                regions_->meshMotionCorrector();

                interfaces_->update();

            }

            regions_->postSolve();

            Info<< "Solved PIMPLE with DNA coupling in "
                << runTime_.cpuTimeIncrement() << " s." << endl;
        }
        else
        {
            //- Solve other partitioned coupled fields
            while (dnaControls_[fldName]->loop())
            {
                assembleAndSolveEqns<fvMatrix, scalar>(fldName);

                assembleAndSolveEqns<fvMatrix, vector>(fldName);

                assembleAndSolveEqns<fvMatrix, tensor>(fldName);

                assembleAndSolveEqns<fvBlockMatrix, vector4>(fldName);

                //assembleAndSolveEqns<symmTensor>(fldName);
            }

            Info<< "Solved "<< fldName << " field with DNA coupling in "
                << runTime_.cpuTimeIncrement() << " s." << endl;
        }
    }


    // Solve region-region coupling (monolithic)
    interfaces_->attach();

    forAll (monolithicCoupledFldNames_, fldI)
    {
        word fldName = monolithicCoupledFldNames_[fldI];

        assembleAndSolveCoupledMatrix<fvMatrix, scalar>
        (
            monolithicCoupledScalarFlds_, fldName
        );

        assembleAndSolveCoupledMatrix<fvMatrix, vector>
        (
            monolithicCoupledVectorFlds_, fldName
        );

        assembleAndSolveCoupledMatrix<fvMatrix, tensor>
        (
            monolithicCoupledTensorFlds_, fldName
        );

//        assembleAndSolveCoupledMatrix<symmTensor>
//        (
//            monolithicCoupledSymmTensorFlds_, fldName
//        );

//        assembleAndSolveCoupledMatrix<fvBlockMatrix, vector4>
//        (
//            monolithicCoupledVector4Flds_, fldName
//        );

        Info<< "Solved "<< fldName << " field with monolithic coupling in "
            << runTime_.cpuTimeIncrement() << " s." << endl;
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
