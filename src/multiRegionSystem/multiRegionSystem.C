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

//! consider HashTable usage 
//template <class T, class Mesh>
//void Foam::fvMeshDistribute::saveBoundaryFields
//(
//    PtrList<FieldField<fvsPatchField, T> >& bflds
//) const
//{
//    typedef const GeometricField<T, fvsPatchField, Mesh> fldType;

//    const HashTable<fldType*> flds
//    (
//        mesh_.objectRegistry::lookupClass<fldType>()
//    );

//    bflds.setSize(flds.size());

//    label i = 0;

//    forAllConstIter (typename HashTable<const fldType*>, flds, iter)
//    {
//        const fldType& fld = *iter();

//        bflds.set(i, fld.boundaryField().clone().ptr());

//        i++;
//    }
//}


template<class T>
void Foam::multiRegionSystem::assembleAndSolveBlockMatrix
(
    PtrList<GeometricField<T, fvPatchField, volMesh> >& flds,
    word fldName
) const
{
    label nEqns = 0;

    forAll (flds, fldI)
    {
//        using fvM = fvMatrix<T>;
        if(flds[fldI].name() == fldName)
        {
            nEqns++;
        }
    }

    Info << "Number of coupled fields : " << nEqns << endl;

    coupledFvScalarMatrix coupledEqns(nEqns);

    forAll (flds, fldI)
    {
        const fvMatrix<T>& eqn = 
            flds[fldI].mesh().objectRegistry::lookupObject<IOReferencer<fvMatrix<T> > >
            (
                fldName + "Eqn"
            )();

        // Add region equation
        Info<< "Adding matrix " 
            << fldName + "Eqn" << " to blockMatrix" 
            << endl;

        // coupledEqns.set(fldI, eqn);

        coupledEqns.set(fldI, &const_cast<fvMatrix<T>& >(eqn)); //! check

        // fvMatrix<T>& eqnRef = const_cast<fvMatrix<T>& >(eqn);
        // coupledEqns.set(fldI, eqnRef);
    }

    coupledEqns.solve(mesh_.solutionDict().solver(fldName + "coupled"));
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

    interfaces_->attach();

    // Get unique list of coupled field names
    hashedWordList fldNames;

    forAllConstIter(fieldsTable, interfaces_->partitionedCoupledFields(), iter)
    {
        // Info << iter() << endl;

        forAll(iter(), fldNameI)
        {
            word fldName = iter()[fldNameI];

            if (!fldNames.contains(fldName))
            {
                fldNames.append(fldName);
            }
        }
    }

    Info << "List of coupled field names : " << fldNames << endl;

    // Assemble coupled fields by type
    PtrList<volScalarField> scalarFlds(0);
    PtrList<volVectorField> vectorFlds(0);
    PtrList<volSphericalTensorField> sphericalTensorFlds(0);
    PtrList<volSymmTensorField> symmTensorFlds(0);
    PtrList<volTensorField> tensorFlds(0);

    assembleCoupledFields<scalar>(scalarFlds, regions(), fldNames);

    // Assemble and solve coupled block-matrix systems
    forAll (fldNames, fldI)
    {
        word fldName = fldNames[fldI];

        assembleAndSolveBlockMatrix<scalar>(scalarFlds, fldName);
    }










////    labelList Ncflds(fldNames.size(), 0);
////    DynamicList<word> Tclfds;

//    const label nVolFieldTypes = 5;
//    const word volFieldTypes[] =
//    {
//        volScalarField::typeName,
//        volVectorField::typeName,
//        volSphericalTensorField::typeName,
//        volSymmTensorField::typeName,
//        volTensorField::typeName
//    };

//    label nScalarFields = 0;
//    label nVectorFlds = 0;
//    label nSphericalTensorFlds = 0;
//    label nSymmTensorFlds = 0;
//    label nTensorFlds = 0;

//    PtrList<volScalarField> scalarFlds(0);
//    PtrList<volVectorField> vectorFlds(0);
//    PtrList<volSphericalTensorField> sphericalTensorFlds(0);
//    PtrList<volSymmTensorField> symmTensorFlds(0);
//    PtrList<volTensorField> tensorFlds(0);

//    forAll (regions(), regI) // go through all regions
//    {
//        // get list of all objects registered to region
//        IOobjectList objects
//        (
//            regions()[regI],
//            "0"
////            regions()[regI].time().timeName()
//        );

//        for (label i=0; i<nVolFieldTypes; i++)
//        {
//            // Search list of field objects for wanted type
//            IOobjectList volTypeObjects = 
//                objects.lookupClass(volFieldTypes[i]);

//            for
//            (
//                IOobjectList::iterator iter = volTypeObjects.begin();
//                iter != volTypeObjects.end();
//                ++iter
//            )
//            {
//                if
//                (
//                    // make sure to pick up only coupled flds of a type
//                    fldNames.contains(iter()->name())
//                 && volFieldTypes[i] == volScalarField::typeName
//                )
//                {
//                    const volScalarField& fld =
//                        regions()[regI].thisDb().lookupObject
//                        <volScalarField>
//                        (
//                            iter()->name()
//                        );

//                    scalarFlds.setSize(scalarFlds.size() + 1);

//                    scalarFlds.set
//                    (
//                        nScalarFields,
//                        fld
//                    );

//                    nScalarFields++;

//                    Info<< "Coupled field type : " << volFieldTypes[i]
//                        << ", name : " << iter()->name()
//                        << endl;
//                }
//            }
//        }
//    }


//    forAll (fldNames, fldI)
//    {
//        word fldName = fldNames[fldI];

////        forAll (regions(), regI) // go through all regions
////        {

//////            for
//////            (
//////                IOobjectList::iterator iter = objects.begin();
//////                iter != objects.end();
//////                ++iter
//////            )
//////            {
//////    //            if (iter()->psi().typeName() == volScalarField::typeName)
//////    //            {
//////    //                Info << "Eqn type : fvScalarMatrix" << endl; 
//////    //            }
//////            }

////        }


//        //-- get number of coupled equations by name
//        label nScalarEqns = 0;
////        IOobjectList eqns(0);

//        forAll (regions(), regI) // go through all regions
//        {
//            if
//            (
//                regions()[regI].foundObject<volScalarField>
//                (
//                    fldNames[fldI]
//                )
//            )
//            {
//                nScalarEqns++;
//            }

//        }

//        //-- assemble block-matrix system
//        PtrList<fvScalarMatrix> scalarEqns(nScalarEqns);

//        label i = 0;

//        forAll (regions(), regI) // go through all regions
//        {
//            if
//            (
//                regions()[regI].foundObject<regIOobject>
//                (
//                    fldNames[fldI] + "Eqn"
//                )
//            )
//            {
//                const fvScalarMatrix& eqn =
//                    regions()[regI].thisDb().lookupObject
//                    <IOReferencer<fvScalarMatrix> >
//                    (
//                        fldNames[fldI] + "Eqn"
//                    )();

//                scalarEqns.set
//                (
//                    i,
//                    eqn
//                );

//                i++;

//                Info<< "Number of coupled matrices : " 
//                    << scalarEqns.size() << endl;

////                const regIOobject& eqn = 
////                    regions()[regI].thisDb().lookupObject<regIOobject>
////                    (
////                        fldNames[fldI] + "Eqn"
////                    );

////                eqns.add(const_cast<regIOobject&>(eqn));

////                Info << "Governing Eqns : " << eqns.toc() << endl;
//            }
//        }

//        if (nScalarEqns =! i)
//        {
//            // Error
//        }


//        // ...



//        for
//        (
//            IOobjectList::iterator iter = eqns.begin();
//            iter != eqns.end();
//            ++iter
//        )
//        {
//            Info << iter()->name() << endl;
////            if (iter()->psi().typeName() == volScalarField::typeName)
////            {
////                Info << "Eqn type : fvScalarMatrix" << endl; 
////            }
//        }
//    }







//    forAll (fldNames, fldI)
//    {
//        word fldName = fldNames[fldI];

//        IOobjectList eqns(0);

//        forAll (regions(), regI) // go through all regions
//        {
//            //! src/foam/fields/ReadFields

//            // get list of all objects registered to region
////            IOobjectList objects
////            (
////                regions()[regI],
////                "0"
//////                regions()[regI].time().timeName()
////            );

////            Info << objects.names() << endl;


//            if
//            (
//                regions()[regI].foundObject<regIOobject>
//                (
//                    fldNames[fldI] + "Eqn"
//                )
//            )
//            {
//                Info << "Found!" << endl;

////                const fvScalarMatrix& eqn =
////                    regions()[regI].thisDb().lookupObject
////                    <IOReferencer<fvScalarMatrix> >
////                    (
////                        fldNames[fldI] + "Eqn"
////                    )();

////                Info << objects.names() << endl;

////                fvScalarMatrix& eqnPtr = const_cast<fvScalarMatrix&>(eqn);
////                forAll (objects, objI)
////                {
////                    if (objects.names()[objI] == fldNames[fldI])
////                    {
////                        const fvScalarMatrix& eqn =
////                            regions()[regI].thisDb().lookupObject
////                            <IOReferencer<fvScalarMatrix> >
////                            (
////                                fldNames[fldI] + "Eqn"
////                            )();

//                        const regIOobject& eqn = 
//                            regions()[regI].thisDb().lookupObject<regIOobject>
//                            (
//                                fldNames[fldI] + "Eqn"
//                            );

//                        eqns.add(const_cast<regIOobject&>(eqn));
////                    }
////                }

//                Info << "eqns : " << eqns.toc() << endl;
//            }

////            for (label i=0; i<nVolFieldTypes; i++)
////            {
////                // Search list of objects for wanted type
////                IOobjectList volTypeObjects = 
////                    objects.lookupClass(volFieldTypes[i]);

////                Info << "field type : " << volTypeObjects.names() << endl;

////                Info << "vol field type : " << volFieldTypes[i] << endl;

//////                const IOobject& eqn = *objects[(fldI + "Eqn")];

//////                const regIOobject& sEqn = 
//////                    regions()[regI].lookupObject<regIOobject>
//////                    ("TEqn");

////                if
////                (
////                    regions()[regI].foundObject<regIOobject>
////                    (
////                        fldNames[fldI] + "Eqn"
////                    )
////                )
////                {
////                    Info << "Found!" << endl;

////                    const fvScalarMatrix& eqn =
////                        regions()[regI].thisDb().lookupObject
////                        <IOReferencer<fvScalarMatrix> >
////                        (
////                            fldNames[fldI] + "Eqn"
////                        )();

////                    //fvScalarMatrix& eqnPtr = const_cast<fvScalarMatrix&>(eqn);
////                }
////            }
//        }
//    }

//    Info << "Number of coupled Fields : " << Ncflds << endl;

    //- assemble matrix systems


//    coupledFvScalarMatrix scalarEqns(NEqns);

//    forAllConstIter(fieldsTable, interfaces_->partitionedCoupledFields(), iter)
//    {
//        forAll(iter(), cfldI)
//        {
//            // Info << iter() << endl;

//            forAll(iter(), fldNameI)
//            {
//                word fldName = iter()[fldNameI];

//                forAll(regions(), regI)
//                {
//                    if (regions()[regI].foundObject<volScalarField>(fldName))
//                    {
//                        Info<< "Found field " << fldName
//                            << " in " << regions()[regI].name()
//                            << " for monolithic coupling"
//                            << endl;

//                        NEqns++;
//                        // scalarEqns.set
//                        // (
//                        //     regI, 
//                        //     regions()[regI].coupledScalarEquations(fldName)
//                        // );
//                    }
//                    // TODO: Sanity check
//                    // else
//                    // {
//                    //     forAll(interfaces(), intI)
//                    //     {
//                    //         FatalError << "Coupling field " << fldName
//                    //             << " not found"
//                    //             << exit(FatalError);
//                    //     }
//                    // }
//                }

//                // scalarEqns.solve
//                // (
//                //     mesh_.solutionDict().solver(fldNameI + "coupled")
//                // );
//            }
//        }
//    }


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
