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
            if (iter().isDict())
            {
                fieldData fd;

                const dictionary& subFieldDict(iter().dict());
                fd.name = iter().keyword();
                fd.maxJumpRes = readScalar(subFieldDict.lookup("maxJumpRes"));
                fd.maxFluxRes = readScalar(subFieldDict.lookup("maxFluxRes"));

                fd.outputJumpResField =
                    subFieldDict.lookupOrDefault<Switch>
                    ("outputJumpResField", false);

                fd.outputFluxResField =
                    subFieldDict.lookupOrDefault<Switch>
                    ("outputFluxResField", false);


                data.append(fd);
            }
        }

        if
        (
            residualDict.found("maxJumpRes")
         && residualDict.found("maxFluxRes")
        )
        {
            fieldData fd;

            fd.name = fldName_;
            fd.maxJumpRes = readScalar(residualDict.lookup("maxJumpRes"));
            fd.maxFluxRes = readScalar(residualDict.lookup("maxFluxRes"));

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

    globalMaxJumpRes_.setSize(dnaResControl_.size());
    globalMaxFluxRes_.setSize(dnaResControl_.size());

    if (debug)
    {
        forAll (dnaResControl_, i)
        {
            const fieldData& fd = dnaResControl_[i];
            Info<< "residualControl[" << i << "]:" << nl
                << "    name         : " << fd.name << nl
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
    scalar& globalMaxJumpRes,
    scalar& globalMaxFluxRes
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    // Check if field of type fieldType and name fldName exists
    // on both sides
    if
    (
        interface.meshA().foundObject<fieldType>(fldName)
     && interface.meshB().foundObject<fieldType>(fldName)
    )
    {
        const fvPatchField<Type>& patchAField =
            interface.patchA().lookupPatchField<fieldType, Type>(fldName);

        const fvPatchField<Type>& patchBField =
            interface.patchB().lookupPatchField<fieldType, Type>(fldName);

        if
        (
            isA<genericRegionCoupledJumpFvPatchField<Type>>(patchAField) &&
            isA<genericRegionCoupledFluxFvPatchField<Type>>(patchBField)
        )
        {
            scalar patchMaxJumpRes =
                refCast<const genericRegionCoupledJumpFvPatchField<Type> >
                (patchAField).normResidual();

            globalMaxJumpRes =
                max
                (
                    globalMaxJumpRes,
                    patchMaxJumpRes
                );

            scalar patchMaxFluxRes =
                refCast<const genericRegionCoupledFluxFvPatchField<Type> >
                (patchBField).normResidual();

            globalMaxFluxRes =
                max
                (
                    globalMaxFluxRes,
                    patchMaxFluxRes
                );

        }
        else if
        (
            isA<genericRegionCoupledJumpFvPatchField<Type>>(patchBField) &&
            isA<genericRegionCoupledFluxFvPatchField<Type>>(patchAField)
        )
        {

            scalar patchMaxJumpRes =
                refCast<const genericRegionCoupledJumpFvPatchField<Type> >
                (patchBField).normResidual();

            globalMaxJumpRes =
                max
                (
                    globalMaxJumpRes,
                    patchMaxJumpRes
                );

            scalar patchMaxFluxRes =
                refCast<const genericRegionCoupledFluxFvPatchField<Type> >
                (patchAField).normResidual();

            globalMaxFluxRes =
                max
                (
                    globalMaxFluxRes,
                    patchMaxFluxRes
                );

        }
        else
        {
            Warning << "Coupled patchFields of field" << fldName
            << " are not derived from genericRegionCoupledJumpFvPatchField"
            << " or genericRegionCoupledFluxFvPatchField." << nl
            << " No DNA interface residual control possible!" <<nl
            << "Setting globalMaxJumpRes and globalMaxFluxRes to maximum value"
            << endl;

            globalMaxJumpRes = GREAT;
            globalMaxFluxRes = GREAT;
        }

        Info<< interface.interfaceName() << " for field " << fldName << ": " << nl
        << "    globalMaxJumpRes: " << globalMaxJumpRes << nl
        << "    globalMaxFluxRes: " << globalMaxFluxRes << endl;
    }
}

void Foam::dnaControl::maxRes
(
    const regionInterfaceType& interface,
    const word& fldName,
    scalar& globalMaxJumpRes,
    scalar& globalMaxFluxRes
)
{
    maxTypeRes<scalar>
        (interface, fldName, globalMaxJumpRes, globalMaxFluxRes);

    maxTypeRes<vector>
        (interface, fldName, globalMaxJumpRes, globalMaxFluxRes);

    // maxTypeRes<sphericalTensor>(interface, fldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeRes<symmTensor>(interface, fldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeRes<tensor>(interface, fldName, globalMaxJumpRes, globalMaxFluxRes);

    // maxTypeRes<vector2>(interface, fldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeRes<vector4>(interface, fldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeRes<vector6>(interface, fldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeRes<vector8>(interface, fldName, globalMaxJumpRes, globalMaxFluxRes);

    // maxTypeRes<sphericalTensor2>(interface, fldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeRes<sphericalTensor4>(interface, fldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeRes<sphericalTensor6>(interface, fldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeRes<sphericalTensor8>(interface, fldName, globalMaxJumpRes, globalMaxFluxRes);

    // maxTypeRes<tensor2>(interface, fldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeRes<tensor4>(interface, fldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeRes<tensor6>(interface, fldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeRes<tensor8>(interface, fldName, globalMaxJumpRes, globalMaxFluxRes);
}

template<class Type>
void Foam::dnaControl::writeResFlds
(
    const regionInterfaceType& interface,
    const word& fldName,
    const Switch outputJumpResField,
    const Switch outputFluxResField,
    bool final
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if
    (
        interface.meshA().foundObject<fieldType>(fldName)
     && interface.meshB().foundObject<fieldType>(fldName)
    )
    {
        const fvPatchField<Type>& patchAField =
            interface.patchA().lookupPatchField<fieldType, Type>(fldName);

        const fvPatchField<Type>& patchBField =
            interface.patchB().lookupPatchField<fieldType, Type>(fldName);

        label patchAID = interface.patchAID();

        if
        (
            isA<genericRegionCoupledJumpFvPatchField<Type>>(patchAField) &&
            isA<genericRegionCoupledFluxFvPatchField<Type>>(patchBField)
        )
        {
            const genericRegionCoupledJumpFvPatchField<Type>& jumpPatchAField =
                refCast<const genericRegionCoupledJumpFvPatchField<Type> >
                (patchAField);

            const genericRegionCoupledFluxFvPatchField<Type>& fluxPatchBField =
                refCast<const genericRegionCoupledFluxFvPatchField<Type> >
                (patchBField);

            if (outputJumpResField)
            {
                if (final)
                {
                    word jumpResFldName = fldName + "finalJumpRes";

                    finalJumpResFlds_[jumpResFldName]->boundaryField()[patchAID] =
                        jumpPatchAField.rawResidual();
                }
                else
                {
                    word jumpResFldName = fldName + "initialJumpRes";

                    initJumpResFlds_[jumpResFldName]->boundaryField()[patchAID] =
                        jumpPatchAField.rawResidual();
                }
            }

            if (outputFluxResField)
            {
                if (final)
                {
                    word fluxResFldName = fldName + "finalFluxRes";

                    tmp<scalarField> resField =
                        jumpPatchAField.interpolateFromNbrField
                        (
                            fluxPatchBField.rawResidual()
                        );

                    finalFluxResFlds_[fluxResFldName]->boundaryField()[patchAID] =
                        resField();
                }
                else
                {
                    word fluxResFldName = fldName + "initialFluxRes";

                    tmp<scalarField> resField =
                        jumpPatchAField.interpolateFromNbrField
                        (
                            fluxPatchBField.rawResidual()
                        );

                    initFluxResFlds_[fluxResFldName]->boundaryField()[patchAID] =
                        resField();
                }
            }

        }
        else if
        (
            isA<genericRegionCoupledJumpFvPatchField<Type>>(patchBField) &&
            isA<genericRegionCoupledFluxFvPatchField<Type>>(patchAField)
        )
        {
            const genericRegionCoupledFluxFvPatchField<Type>& fluxPatchAField =
                refCast<const genericRegionCoupledFluxFvPatchField<Type> >
                (patchAField);

            const genericRegionCoupledJumpFvPatchField<Type>& jumpPatchBField =
                refCast<const genericRegionCoupledJumpFvPatchField<Type> >
                (patchBField);

            if (outputFluxResField)
            {
                if (final)
                {
                    word fluxResFldName = fldName + "finalFluxRes";

                    finalFluxResFlds_[fluxResFldName]->boundaryField()[patchAID] =
                        fluxPatchAField.rawResidual();
                }
                else
                {
                    word fluxResFldName = fldName + "initialFluxRes";

                    initFluxResFlds_[fluxResFldName]->boundaryField()[patchAID] =
                        fluxPatchAField.rawResidual();
                }
            }

            if (outputJumpResField)
            {
                if (final)
                {
                    word jumpResFldName = fldName + "finalJumpRes";

                    tmp<scalarField> resField =
                        fluxPatchAField.interpolateFromNbrField
                        (
                            jumpPatchBField.rawResidual()
                        );

                    finalJumpResFlds_[jumpResFldName]->boundaryField()[patchAID] =
                        resField();
                }
                else
                {
                    word jumpResFldName = fldName + "initialJumpRes";

                    tmp<scalarField> resField =
                        fluxPatchAField.interpolateFromNbrField
                        (
                            jumpPatchBField.rawResidual()
                        );

                    initJumpResFlds_[jumpResFldName]->boundaryField()[patchAID] =
                        resField();
                }
            }
        }
    }
}

void Foam::dnaControl::writeResFlds
(
    const regionInterfaceType& interface,
    const word& fldName,
    const Switch outputJumpResField,
    const Switch outputFluxResField,
    bool final
)
{
    writeResFlds<scalar>
        (interface, fldName, outputJumpResField, outputFluxResField, final);

    writeResFlds<vector>
        (interface, fldName, outputJumpResField, outputFluxResField, final);

    // writeResFlds<sphericalTensor>(interface, fldName, outputJumpResField, outputFluxResField, final);
    // writeResFlds<symmTensor>(interface, fldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFlds<tensor>(interface, fldName, outputJumpResFieldoutputFluxResField, final);

    // writeResFlds<vector2>(interface, fldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFlds<vector4>(interface, fldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFlds<vector6>(interface, fldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFlds<vector8>(interface, fldName, outputJumpResFieldoutputFluxResField, final);

    // writeResFlds<sphericalTensor2>(interface, fldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFlds<sphericalTensor4>(interface, fldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFlds<sphericalTensor6>(interface, fldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFlds<sphericalTensor8>(interface, fldName, outputJumpResFieldoutputFluxResField, final);

    // writeResFlds<tensor2>(interface, fldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFlds<tensor4>(interface, fldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFlds<tensor6>(interface, fldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFlds<tensor8>(interface, fldName, outputJumpResFieldoutputFluxResField, final);
}

void Foam::dnaControl::outputMaxResInfo
(
    bool criteriaSatisfied,
    label corr
)
{
    if ((!dnaResControl_.empty()))
    {
        // output initial and final residual
        if(corr == 2)
        {
            forAll (dnaResControl_, ctrlI)
            {
                word fldName = dnaResControl_[ctrlI].name;

                Info<< "    "
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
                    word subFieldName = dnaResControl_[ctrlI].name;

                    Switch outputJumpResField =
                        dnaResControl_[ctrlI].outputJumpResField;

                    Switch outputFluxResField =
                        dnaResControl_[ctrlI].outputFluxResField;

                    writeResFlds
                    (
                        interfaces_[intI],
                        subFieldName,
                        outputJumpResField,
                        outputFluxResField,
                        false
                    );
                }
            }
        }

        if(criteriaSatisfied || corr > maxCoupleIter_)
        {
            forAll (dnaResControl_, ctrlI)
            {
                word fldName = dnaResControl_[ctrlI].name;

                Info<< "    "
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
                    word subFieldName = dnaResControl_[ctrlI].name;

                    Switch outputJumpResField =
                        dnaResControl_[ctrlI].outputJumpResField;

                    Switch outputFluxResField =
                        dnaResControl_[ctrlI].outputFluxResField;

                    writeResFlds
                    (
                        interfaces_[intI],
                        subFieldName,
                        outputJumpResField,
                        outputFluxResField,
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
    globalMaxJumpRes_ = 0;
    globalMaxFluxRes_ = 0;

    forAll (interfaces_, intI)
    {
        forAll (dnaResControl_, ctrlI)
        {
            word subFieldName = dnaResControl_[ctrlI].name;

            maxRes
            (
                interfaces_[intI],
                subFieldName,
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
         && (globalMaxJumpRes_[ctrlI] <= dnaResControl_[ctrlI].maxJumpRes);

        criteriaSatisfied =
            criteriaSatisfied
         && (globalMaxFluxRes_[ctrlI] <= dnaResControl_[ctrlI].maxFluxRes);
    }



    return criteriaSatisfied;

    // return false;
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
                // initialJumpResFields
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

                // finalJumpResFields
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
    globalMaxJumpRes_(),
    globalMaxFluxRes_(),
    initJumpResFlds_(),
    initFluxResFlds_(),
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
				Info<< "    field " << dnaResControl_[i].name << token::TAB
					<< ": maxJumpRes " << dnaResControl_[i].maxJumpRes
					<< ", maxFluxRes " << dnaResControl_[i].maxFluxRes
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

        corr_ = 0;

        return false;
    }

    Info<< nl
    << "DNA iteration: " << corr_
    << endl;

    return true;
}
