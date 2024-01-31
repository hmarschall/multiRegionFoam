/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dnaControl.H"
#include "fieldTypes.H"
#include "VectorNFieldTypes.H"
#include "scalar.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dnaControl, 0);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::dnaControl::read()
{
    const dictionary& dnaFieldDict = this->dict().subDict(fldName_);

    maxCoupleIter_ = dnaFieldDict.lookupOrDefault<label>("maxCoupleIter", 1);

    if (dnaFieldDict.found("residualControl"))
    {
        const dictionary& residualDict = dnaFieldDict.subDict("residualControl");

        DynamicList<fieldData> data(0);

        forAllConstIter(dictionary, residualDict, iter)
        {
            // If multiple fields are coupled through the interface at once
            // e.g. p-U coupling
            if (iter().isDict())
            {
                fieldData fd;

                const dictionary& subFieldDict(iter().dict());
                fd.name = iter().keyword();

                if
                (
                    subFieldDict.found("maxIntRes")
                 && subFieldDict.found("maxJumpRes")
                 && subFieldDict.found("maxFluxRes")
                )
                {
                    FatalErrorIn
                    (
                        "void Foam::dnaControl::read()"
                    )
                    << "Field you can only specify either a single maxIntRes "
                    << "residual or both the maxJumpRes and maxFluxRes"
                    << abort(FatalError);
                }
                else if
                (
                    subFieldDict.found("maxIntRes")
                )
                {
                    // If only the single interface residual is specified
                    // set the jump and flux residual to be satisfied by
                    // default
                    fd.maxIntRes = readScalar(subFieldDict.lookup("maxIntRes"));
                    fd.maxJumpRes = VGREAT;
                    fd.maxFluxRes = VGREAT;

                    fd.outputIntResField =
                        subFieldDict.lookupOrDefault<Switch>
                        ("outputIntResField", false);
                    fd.outputJumpResField = false;
                    fd.outputFluxResField = false;
                }
                else if
                (
                    subFieldDict.found("maxJumpRes")
                 && subFieldDict.found("maxFluxRes")
                )
                {
                    // If only the jump and flux residual are specified
                    // set the single residual to be satisfied by
                    // default
                    fd.maxIntRes = VGREAT;
                    fd.maxJumpRes = readScalar(subFieldDict.lookup("maxJumpRes"));
                    fd.maxFluxRes = readScalar(subFieldDict.lookup("maxFluxRes"));

                    fd.outputIntResField = false;
                    fd.outputJumpResField =
                        subFieldDict.lookupOrDefault<Switch>
                        ("outputJumpResField", false);
                    fd.outputFluxResField =
                        subFieldDict.lookupOrDefault<Switch>
                        ("outputFluxResField", false);
                }

                data.append(fd);
            }
        }

        if
        (
            residualDict.found("maxIntRes")
         && residualDict.found("maxJumpRes")
         && residualDict.found("maxFluxRes")
        )
        {
            FatalErrorIn
            (
                "void Foam::dnaControl::read()"
            )
            << "Field you can only specify either a single maxIntRes "
            << "residual or both the maxJumpRes and maxFluxRes"
            << abort(FatalError);
        }
        else if
        (
            residualDict.found("maxIntRes")
        )
        {
            fieldData fd;
            fd.name = fldName_;
            // If only the single interface residual is specified
            // set the jump and flux residual to be satisfied by
            // default
            fd.maxIntRes = readScalar(residualDict.lookup("maxIntRes"));
            fd.maxJumpRes = VGREAT;
            fd.maxFluxRes = VGREAT;

            fd.outputIntResField =
                residualDict.lookupOrDefault<Switch>
                ("outputIntResField", false);
            fd.outputJumpResField = false;
            fd.outputFluxResField = false;

            data.append(fd);
        }
        else if
        (
            residualDict.found("maxJumpRes")
         && residualDict.found("maxFluxRes")
        )
        {
            fieldData fd;
            fd.name = fldName_;
            // If only the jump and flux residual are specified
            // set the single residual to be satisfied by
            // default
            fd.maxIntRes = VGREAT;
            fd.maxJumpRes = readScalar(residualDict.lookup("maxJumpRes"));
            fd.maxFluxRes = readScalar(residualDict.lookup("maxFluxRes"));

            fd.outputIntResField = false;
            fd.outputJumpResField =
                residualDict.lookupOrDefault<Switch>
                ("outputJumpResField", false);
            fd.outputFluxResField =
                residualDict.lookupOrDefault<Switch>
                ("outputFluxResField", false);

            data.append(fd);
        }

        dnaResControl_.transfer(data);
    }

    globalMaxIntRes_.setSize(dnaResControl_.size());
    globalMaxJumpRes_.setSize(dnaResControl_.size());
    globalMaxFluxRes_.setSize(dnaResControl_.size());

    if (debug)
    {
        forAll (dnaResControl_, i)
        {
            const fieldData& fd = dnaResControl_[i];
            Info<< "residualControl[" << i << "]:" << nl
                << "    name         : " << fd.name << nl
                << "    maxIntRes   : " << fd.maxIntRes << nl
                << "    maxJumpRes   : " << fd.maxJumpRes << nl
                << "    maxFluxRes   : " << fd.maxFluxRes << endl;
        }
    }
}

template<class Type>
void Foam::dnaControl::maxTypeRes
(
    const regionInterfaceType& interface,
    const word& fldName,
    scalar& globalMaxIntRes,
    scalar& globalMaxJumpRes,
    scalar& globalMaxFluxRes
)
{
    scalar patchMaxIntRes = interface.normIntRes(fldName);

    globalMaxIntRes =
        max
        (
            globalMaxIntRes,
            patchMaxIntRes
        );

    if (interface.couplesField<Type>(fldName))
    {
        scalar patchMaxJumpRes = interface.normJumpRes<Type>(fldName);
        scalar patchMaxFluxRes = interface.normFluxRes<Type>(fldName);

        globalMaxJumpRes =
            max
            (
                globalMaxJumpRes,
                patchMaxJumpRes
            );

        globalMaxFluxRes =
            max
            (
                globalMaxFluxRes,
                patchMaxFluxRes
            );
    }

    Info<< interface.interfaceName() << " for field " << fldName << ": " << nl
        << "    globalMaxIntRes: " << globalMaxIntRes << nl
        << "    globalMaxJumpRes: " << globalMaxJumpRes << nl
        << "    globalMaxFluxRes: " << globalMaxFluxRes << endl;
}

void Foam::dnaControl::maxRes
(
    const regionInterfaceType& interface,
    const word& fldName,
    scalar& globalMaxIntRes,
    scalar& globalMaxJumpRes,
    scalar& globalMaxFluxRes
)
{
    maxTypeRes<scalar>
    (
        interface, fldName,globalMaxIntRes, globalMaxJumpRes, globalMaxFluxRes
    );
    maxTypeRes<vector>
    (
        interface, fldName,globalMaxIntRes, globalMaxJumpRes, globalMaxFluxRes
    );

    // maxTypeRes<sphericalTensor>
    // (
    //     interface, fldName,globalMaxIntRes, globalMaxJumpRes, globalMaxFluxRes
    // );
    // maxTypeRes<symmTensor>
    // (
    //     interface, fldName,globalMaxIntRes, globalMaxJumpRes, globalMaxFluxRes
    // );
    // maxTypeRes<tensor>
    // (
    //     interface, fldName,globalMaxIntRes, globalMaxJumpRes, globalMaxFluxRes
    // );

    // maxTypeRes<vector2>
    // (
    //     interface, fldName,globalMaxIntRes, globalMaxJumpRes, globalMaxFluxRes
    // );
    // maxTypeRes<vector4>
    // (
    //     interface, fldName,globalMaxIntRes, globalMaxJumpRes, globalMaxFluxRes
    // );
    // maxTypeRes<vector6>
    // (
    //     interface, fldName,globalMaxIntRes, globalMaxJumpRes, globalMaxFluxRes
    // );
    // maxTypeRes<vector8>
    // (
    //     interface, fldName,globalMaxIntRes, globalMaxJumpRes, globalMaxFluxRes
    // );

    // maxTypeRes<sphericalTensor2>
    // (
    //     interface, fldName,globalMaxIntRes, globalMaxJumpRes, globalMaxFluxRes
    // );
    // maxTypeRes<sphericalTensor4>
    // (
    //     interface, fldName,globalMaxIntRes, globalMaxJumpRes, globalMaxFluxRes
    // );
    // maxTypeRes<sphericalTensor6>
    // (
    //     interface, fldName,globalMaxIntRes, globalMaxJumpRes, globalMaxFluxRes
    // );
    // maxTypeRes<sphericalTensor8>
    // (
    //     interface, fldName,globalMaxIntRes, globalMaxJumpRes, globalMaxFluxRes
    // );

    // maxTypeRes<tensor2>
    // (
    //     interface, fldName,globalMaxIntRes, globalMaxJumpRes, globalMaxFluxRes
    // );
    // maxTypeRes<tensor4>
    // (
    //     interface, fldName,globalMaxIntRes, globalMaxJumpRes, globalMaxFluxRes
    // );
    // maxTypeRes<tensor6>
    // (
    //     interface, fldName,globalMaxIntRes, globalMaxJumpRes, globalMaxFluxRes
    // );
    // maxTypeRes<tensor8>
    // (
    //     interface, fldName,globalMaxIntRes, globalMaxJumpRes, globalMaxFluxRes
    // );
}

template<class Type>
void Foam::dnaControl::writeResFlds
(
    const regionInterfaceType& interface,
    const word& fldName,
    const Switch outputIntResField,
    const Switch outputJumpResField,
    const Switch outputFluxResField,
    bool final
)
{
    label patchAID = interface.patchAID();

    if (outputIntResField)
    {
        if (final)
        {
            word intResFldName = fldName + "finalIntRes";

            finalIntResFlds_[intResFldName]->boundaryField()[patchAID] =
                interface.rawIntRes(fldName);
        }
        else
        {
            word intResFldName = fldName + "initialIntRes";

            initIntResFlds_[intResFldName]->boundaryField()[patchAID] =
                interface.rawIntRes(fldName);
        }
    }

    if (interface.couplesField<Type>(fldName))
    {
        if (outputJumpResField)
        {
            if (final)
            {
                word jumpResFldName = fldName + "finalJumpRes";

                finalJumpResFlds_[jumpResFldName]->boundaryField()[patchAID] =
                    interface.rawJumpRes<Type>(fldName);
            }
            else
            {
                word jumpResFldName = fldName + "initialJumpRes";

                initJumpResFlds_[jumpResFldName]->boundaryField()[patchAID] =
                    interface.rawJumpRes<Type>(fldName);
            }
        }

        if (outputFluxResField)
        {
            if (final)
            {
                word fluxResFldName = fldName + "finalFluxRes";

                finalFluxResFlds_[fluxResFldName]->boundaryField()[patchAID] =
                    interface.rawFluxRes<Type>(fldName);
            }
            else
            {
                word fluxResFldName = fldName + "initialFluxRes";

                initFluxResFlds_[fluxResFldName]->boundaryField()[patchAID] =
                    interface.rawFluxRes<Type>(fldName);
            }
        }
    }
}

void Foam::dnaControl::writeResFlds
(
    const regionInterfaceType& interface,
    const word& fldName,
    const Switch outputIntResField,
    const Switch outputJumpResField,
    const Switch outputFluxResField,
    bool final
)
{
    writeResFlds<scalar>
    (
        interface, fldName, outputIntResField, outputJumpResField, outputFluxResField, final
    );

    writeResFlds<vector>
    (
        interface, fldName, outputIntResField, outputJumpResField, outputFluxResField, final
    );

    // writeResFlds<sphericalTensor>
    // (
    //     interface, fldName, outputIntResField, outputJumpResField, outputFluxResField, final
    // );
    // writeResFlds<symmTensor>
    // (
    //     interface, fldName, outputIntResField, outputJumpResField, outputFluxResField, final
    // );
    // writeResFlds<tensor>
    // (
    //     interface, fldName, outputIntResField, outputJumpResField, outputFluxResField, final
    // );

    // writeResFlds<vector2>
    // (
    //     interface, fldName, outputIntResField, outputJumpResField, outputFluxResField, final
    // );
    // writeResFlds<vector4>
    // (
    //     interface, fldName, outputIntResField, outputJumpResField, outputFluxResField, final
    // );
    // writeResFlds<vector6>
    // (
    //     interface, fldName, outputIntResField, outputJumpResField, outputFluxResField, final
    // );
    // writeResFlds<vector8>
    // (
    //     interface, fldName, outputIntResField, outputJumpResField, outputFluxResField, final
    // );

    // writeResFlds<sphericalTensor2>
    // (
    //     interface, fldName, outputIntResField, outputJumpResField, outputFluxResField, final
    // );
    // writeResFlds<sphericalTensor4>
    // (
    //     interface, fldName, outputIntResField, outputJumpResField, outputFluxResField, final
    // );
    // writeResFlds<sphericalTensor6>
    // (
    //     interface, fldName, outputIntResField, outputJumpResField, outputFluxResField, final
    // );
    // writeResFlds<sphericalTensor8>
    // (
    //     interface, fldName, outputIntResField, outputJumpResField, outputFluxResField, final
    // );

    // writeResFlds<tensor2>
    // (
    //     interface, fldName, outputIntResField, outputJumpResField, outputFluxResField, final
    // );
    // writeResFlds<tensor4>
    // (
    //     interface, fldName, outputIntResField, outputJumpResField, outputFluxResField, final
    // );
    // writeResFlds<tensor6>
    // (
    //     interface, fldName, outputIntResField, outputJumpResField, outputFluxResField, final
    // );
    // writeResFlds<tensor8>
    // (
    //     interface, fldName, outputIntResField, outputJumpResField, outputFluxResField, final
    // );
}

void Foam::dnaControl::outputMaxResInfo
(
    bool criteriaSatisfied,
    label corr
)
{
    if ((!dnaResControl_.empty()))
    {
        // output initial residual
        if(corr == 2)
        {
            forAll (dnaResControl_, ctrlI)
            {
                word fldName = dnaResControl_[ctrlI].name;

                Info<< "\n    "
                    << fldName
                    << " initialMaxIntRes: "
                    << globalMaxIntRes_[ctrlI]
                    << nl
                    << "    "
                    << fldName
                    << " initialMaxJumpRes: "
                    << globalMaxJumpRes_[ctrlI]
                    << nl
                    << "    "
                    << fldName
                    << " initialMaxFluxRes: "
                    << globalMaxFluxRes_[ctrlI]
                    << endl;
            }

            forAll (interfaces_, intI)
            {
                forAll (dnaResControl_, ctrlI)
                {
                    writeResFlds
                    (
                        interfaces_[intI],
                        dnaResControl_[ctrlI].name,
                        dnaResControl_[ctrlI].outputIntResField,
                        dnaResControl_[ctrlI].outputJumpResField,
                        dnaResControl_[ctrlI].outputFluxResField,
                        false
                    );
                }
            }
        }

        // output final residual
        if(criteriaSatisfied || corr > maxCoupleIter_)
        {
            forAll (dnaResControl_, ctrlI)
            {
                word fldName = dnaResControl_[ctrlI].name;

                Info<< "\n    "
                    << fldName
                    << " finalMaxIntRes: "
                    << globalMaxIntRes_[ctrlI]
                    << nl
                    << "    "
                    << fldName
                    << " finalMaxJumpRes: "
                    << globalMaxJumpRes_[ctrlI]
                    << nl
                    << "    "
                    << fldName
                    << " finalMaxFluxRes: "
                    << globalMaxFluxRes_[ctrlI]
                    << endl;
            }

            forAll (interfaces_, intI)
            {
                forAll (dnaResControl_, ctrlI)
                {
                    writeResFlds
                    (
                        interfaces_[intI],
                        dnaResControl_[ctrlI].name,
                        dnaResControl_[ctrlI].outputIntResField,
                        dnaResControl_[ctrlI].outputJumpResField,
                        dnaResControl_[ctrlI].outputFluxResField,
                        true
                    );
                }
            }
        }
    }
}

bool Foam::dnaControl::criteriaSatisfied()
{
    if (dnaResControl_.empty())
    {
        return false;
    }

    // no checks on first iteration - nothing has been calculated yet
    if (corr_ == 1)
    {
        return false;
    }

    //- Reset global max residuals to zero in each coupling iteration
    globalMaxIntRes_ = 0.0;
    globalMaxJumpRes_ = 0.0;
    globalMaxFluxRes_ = 0.0;

    forAll (interfaces_, intI)
    {
        forAll (dnaResControl_, ctrlI)
        {
            maxRes
            (
                interfaces_[intI],
                dnaResControl_[ctrlI].name,
                globalMaxIntRes_[ctrlI],
                globalMaxJumpRes_[ctrlI],
                globalMaxFluxRes_[ctrlI]
            );
        }
    }

    bool criteriaSatisfied = true;

    forAll (dnaResControl_, ctrlI)
    {
        criteriaSatisfied =
            criteriaSatisfied
         && (globalMaxIntRes_[ctrlI] <= dnaResControl_[ctrlI].maxIntRes);

        criteriaSatisfied =
            criteriaSatisfied
         && (globalMaxJumpRes_[ctrlI] <= dnaResControl_[ctrlI].maxJumpRes);

        criteriaSatisfied =
            criteriaSatisfied
         && (globalMaxFluxRes_[ctrlI] <= dnaResControl_[ctrlI].maxFluxRes);
    }

    return criteriaSatisfied;
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dnaControl::createResFlds()
{
    forAll (interfaces_, intI)
    {
        const regionInterfaceType& interface = interfaces_[intI];

        forAll (dnaResControl_, ctrlI)
        {
            fieldData& control = dnaResControl_[ctrlI];

            if (control.outputIntResField)
            {
                // initialIntResFields
                if
                (
                    interface.meshA().foundObject<volScalarField>
                    (control.name + "initialIntRes")
                )
                {
                    initIntResFlds_.set
                    (
                        control.name + "initialIntRes",
                        const_cast<volScalarField*>
                        (
                            &interface.meshA().lookupObject<volScalarField>
                            (control.name + "initialIntRes")
                        )
                    );
                }
                else
                {
                    Info<< "DNA-Control: Creating field "
                        << control.name << "initialIntRes"
                        << endl;

                    initIntResFlds_.set
                    (
                        control.name + "initialIntRes",
                        new volScalarField
                        (
                            IOobject
                            (
                                control.name + "initialIntRes",
                                interface.meshA().time().timeName(),
                                interface.meshA(),
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                            ),
                            interface.meshA(),
                            dimensionedScalar("zero", dimless, 0.)
                        )
                    );
                }

                // finalIntResFields
                if
                (
                    interface.meshA().foundObject<volScalarField>
                    (control.name + "finalIntRes")
                )
                {
                    finalIntResFlds_.set
                    (
                        control.name + "finalIntRes",
                        const_cast<volScalarField*>
                        (
                            &interface.meshA().lookupObject<volScalarField>
                            (control.name + "finalIntRes")
                        )
                    );
                }
                else
                {
                    Info<< "DNA-Control: Creating field "
                        << control.name << "finalIntRes"
                        << endl;

                    finalIntResFlds_.set
                    (
                        control.name + "finalIntRes",
                        new volScalarField
                        (
                            IOobject
                            (
                                control.name + "finalIntRes",
                                interface.meshA().time().timeName(),
                                interface.meshA(),
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                            ),
                            interface.meshA(),
                            dimensionedScalar("zero", dimless, 0.)
                        )
                    );
                }
            }

            if (control.outputJumpResField)
            {
                // initialJumpResFields
                if
                (
                    interface.meshA().foundObject<volScalarField>
                    (control.name + "initialJumpRes")
                )
                {
                    initJumpResFlds_.set
                    (
                        control.name + "initialJumpRes",
                        const_cast<volScalarField*>
                        (
                            &interface.meshA().lookupObject<volScalarField>
                            (control.name + "initialJumpRes")
                        )
                    );
                }
                else
                {
                    Info<< "DNA-Control: Creating field "
                        << control.name << "initialJumpRes"
                        << endl;

                    initJumpResFlds_.set
                    (
                        control.name + "initialJumpRes",
                        new volScalarField
                        (
                            IOobject
                            (
                                control.name + "initialJumpRes",
                                interface.meshA().time().timeName(),
                                interface.meshA(),
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                            ),
                            interface.meshA(),
                            dimensionedScalar("zero", dimless, 0.)
                        )
                    );
                }

                // finalJumpResFields
                if
                (
                    interface.meshA().foundObject<volScalarField>
                    (control.name + "finalJumpRes")
                )
                {
                    finalJumpResFlds_.set
                    (
                        control.name + "finalJumpRes",
                        const_cast<volScalarField*>
                        (
                            &interface.meshA().lookupObject<volScalarField>
                            (control.name + "finalJumpRes")
                        )
                    );
                }
                else
                {
                    Info<< "DNA-Control: Creating field "
                        << control.name << "finalJumpRes"
                        << endl;

                    finalJumpResFlds_.set
                    (
                        control.name + "finalJumpRes",
                        new volScalarField
                        (
                            IOobject
                            (
                                control.name + "finalJumpRes",
                                interface.meshA().time().timeName(),
                                interface.meshA(),
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                            ),
                            interface.meshA(),
                            dimensionedScalar("zero", dimless, 0.)
                        )
                    );
                }
            }

            if (control.outputFluxResField)
            {
                // initialFluxResFields
                if
                (
                    interface.meshA().foundObject<volScalarField>
                    (control.name + "initialFluxRes")
                )
                {
                    initFluxResFlds_.set
                    (
                        control.name + "initialFluxRes",
                        const_cast<volScalarField*>
                        (
                            &interface.meshA().lookupObject<volScalarField>
                            (control.name + "initialFluxRes")
                        )
                    );
                }
                else
                {
                    Info<< "DNA-Control: Creating field "
                        << control.name
                        << "initialFluxRes"
                        << endl;

                    initFluxResFlds_.set
                    (
                        control.name + "initialFluxRes",
                        new volScalarField
                        (
                            IOobject
                            (
                                control.name + "initialFluxRes",
                                interface.meshA().time().timeName(),
                                interface.meshA(),
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                            ),
                            interface.meshA(),
                            dimensionedScalar("zero", dimless, 0.)
                        )
                    );
                }

                // finalFluxResFields
                if
                (
                    interface.meshA().foundObject<volScalarField>
                    (control.name + "finalFluxRes")
                )
                {
                    finalFluxResFlds_.set
                    (
                        control.name + "finalFluxRes",
                        const_cast<volScalarField*>
                        (
                            &interface.meshA().lookupObject<volScalarField>
                            (control.name + "finalFluxRes")
                        )
                    );
                }
                else
                {
                    Info<< "DNA-Control: Creating field "
                        << control.name
                        << "finalFluxRes"
                        << endl;

                    finalFluxResFlds_.set
                    (
                        control.name + "finalFluxRes",
                        new volScalarField
                        (
                            IOobject
                            (
                                control.name + "finalFluxRes",
                                interface.meshA().time().timeName(),
                                interface.meshA(),
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                            ),
                            interface.meshA(),
                            dimensionedScalar("zero", dimless, 0.)
                        )
                    );
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dnaControl::dnaControl
(
    const Time& runTime,
    const word& fldName,
    const regionInterfaceTypeList& interfaces
)
:
    IOobject
    (
        "dnaControl",
        runTime.timeName(),
        runTime
    ),
    fldName_(fldName),
    multiRegionProperties_
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
    interfaces_(interfaces),
    dnaResControl_(),
    globalMaxIntRes_(),
    globalMaxJumpRes_(),
    globalMaxFluxRes_(),
    initIntResFlds_(),
    initJumpResFlds_(),
    initFluxResFlds_(),
    finalIntResFlds_(),
    finalJumpResFlds_(),
    finalFluxResFlds_(),
    maxCoupleIter_(0),
    corr_(0)
{
    read();

    createResFlds();

    Info<< nl;
    if (dnaResControl_.empty())
    {
            Info<< "DNA-Control: no residual control data found. "
                << "Calculations will use "
                << maxCoupleIter_
                << " outer coupling loops"
                << nl << endl;
    }
    else
    {
            Info<< "DNA-Control: max iterations = "
                << maxCoupleIter_
                << endl;

            forAll (dnaResControl_, i)
            {
                Info<< "    field: " << dnaResControl_[i].name << token::TAB
                    << ", maxIntRes: " << dnaResControl_[i].maxIntRes
                    << ", maxJumpRes: " << dnaResControl_[i].maxJumpRes
                    << ", maxFluxRes: " << dnaResControl_[i].maxFluxRes
                    << nl;
            }
            Info<< endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dnaControl::~dnaControl()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::dnaControl::loop()
{
    read();

    corr_++;

    // check if interfaces are coupled yet
    bool interfacesCoupled = false;
    forAll (interfaces_, intI)
    {
        interfacesCoupled = interfacesCoupled || interfaces_[intI].coupled();
    }
    // if all interfaces are not coupled yet -> end the DNA loop
    if
    (
        (!interfacesCoupled) && (corr_>1)
    )
    {
        Info<< nl
            << "Interfaces are not coupled yet" << nl
            << "No DNA coupling is applied"
            << endl;

        corr_ = 0;

        return false;
    }

    bool satisfied = criteriaSatisfied();

    outputMaxResInfo(satisfied, corr_);

    if (satisfied)
    {
        Info<< "DNA: converged in "
            << corr_ - 1
            << " iterations"
            << endl;

        corr_ = 0;

        return false;
    }

    if (corr_ > maxCoupleIter_)
    {
        if ((!dnaResControl_.empty()))
        {
            Info<< "DNA: not converged within "
                << maxCoupleIter_
                << " iterations"
                << endl;
        }
        else
        {
            Info<< "DNA: reached max number of "
                << maxCoupleIter_
                << " iterations"
                << endl;
        }

        corr_ = 0;

        return false;
    }

    Info<< nl
    << "DNA iteration: " << corr_
    << endl;

    return true;
}