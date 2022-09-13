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
    const dictionary& dnaFieldDict = this->dict().subDict(fieldName_);

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

                fd.outputJumpResField = subFieldDict.lookupOrDefault<Switch>("outputJumpResField", false);
                fd.outputFluxResField = subFieldDict.lookupOrDefault<Switch>("outputFluxResField", false);


                data.append(fd);
            }
        }

        if (residualDict.found("maxJumpRes") && residualDict.found("maxFluxRes"))
        {
            fieldData fd;

            fd.name = fieldName_;
            fd.maxJumpRes = readScalar(residualDict.lookup("maxJumpRes"));
            fd.maxFluxRes = readScalar(residualDict.lookup("maxFluxRes"));

            fd.outputJumpResField = residualDict.lookupOrDefault<Switch>("outputJumpResField", false);
            fd.outputFluxResField = residualDict.lookupOrDefault<Switch>("outputFluxResField", false);

            data.append(fd);
        }

        dnaResidualControl_.transfer(data);
    }

    if (debug)
    {
        forAll(dnaResidualControl_, i)
        {
            const fieldData& fd = dnaResidualControl_[i];
            Info<< "residualControl[" << i << "]:" << nl
                << "    name         : " << fd.name << nl
                << "    maxJumpRes   : " << fd.maxJumpRes << nl
                << "    maxFluxRes   : " << fd.maxFluxRes << endl;
        }
    }
}

template<class Type>
void Foam::dnaControl::maxTypeResidual
(
    const regionInterface& interface,
    const word& fieldName,
    scalar& globalMaxJumpRes,
    scalar& globalMaxFluxRes
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;
    
    // Check if field of type fieldType and name fieldName exists
    // - since the field is a coupled field it is assumed that 
    //   if it exists on meshA it also exists on meshB 
    if (interface.meshA().foundObject<fieldType>(fieldName))
    {
        const fvPatchField<Type>& patchAField =
            interface.patchA().lookupPatchField<fieldType, Type>(fieldName);

        const fvPatchField<Type>& patchBField =
            interface.patchB().lookupPatchField<fieldType, Type>(fieldName);

        if 
        (
            isA<genericRegionCoupledJumpFvPatchField<Type>>(patchAField) &&
            isA<genericRegionCoupledFluxFvPatchField<Type>>(patchBField)
        )
        {
            scalar patchMaxJumpRes = refCast<const genericRegionCoupledJumpFvPatchField<Type>>(patchAField).ofNormResidual();
            globalMaxJumpRes = 
                max
                (
                    globalMaxJumpRes,
                    patchMaxJumpRes
                );
            
            scalar patchMaxFluxRes = refCast<const genericRegionCoupledFluxFvPatchField<Type>>(patchBField).ofNormResidual();
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
            
            scalar patchMaxJumpRes = refCast<const genericRegionCoupledJumpFvPatchField<Type>>(patchBField).ofNormResidual();
            globalMaxJumpRes = 
                max
                (
                    globalMaxJumpRes,
                    patchMaxJumpRes
                );
            
            scalar patchMaxFluxRes = refCast<const genericRegionCoupledFluxFvPatchField<Type>>(patchAField).ofNormResidual();
            globalMaxFluxRes = 
                max
                (
                    globalMaxFluxRes,
                    patchMaxFluxRes
                );

        }
        else
        {
            Warning << "Coupled patchFields of field" << fieldName 
            << " are not derived from genericRegionCoupledJumpFvPatchField" 
            << " or genericRegionCoupledFluxFvPatchField." << nl
            << " No DNA interface residual control possible!" <<nl
            << "Setting globalMaxJumpRes and globalMaxFluxRes to maximum value"
            << endl;

            globalMaxJumpRes = GREAT;
            globalMaxFluxRes = GREAT;
        }

        Info<< interface.name() << " for field " << fieldName << ": " << nl
        << "    globalMaxJumpRes: " << globalMaxJumpRes << nl 
        << "    globalMaxFluxRes: " << globalMaxFluxRes << endl;
    }
}

void Foam::dnaControl::maxResidual
(
    const regionInterface& interface,
    const word& fieldName,
    scalar& globalMaxJumpRes,
    scalar& globalMaxFluxRes
)
{
    maxTypeResidual<scalar>(interface, fieldName, globalMaxJumpRes, globalMaxFluxRes);
    maxTypeResidual<vector>(interface, fieldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeResidual<sphericalTensor>(interface, fieldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeResidual<symmTensor>(interface, fieldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeResidual<tensor>(interface, fieldName, globalMaxJumpRes, globalMaxFluxRes);

    // maxTypeResidual<vector2>(interface, fieldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeResidual<vector4>(interface, fieldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeResidual<vector6>(interface, fieldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeResidual<vector8>(interface, fieldName, globalMaxJumpRes, globalMaxFluxRes);

    // maxTypeResidual<sphericalTensor2>(interface, fieldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeResidual<sphericalTensor4>(interface, fieldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeResidual<sphericalTensor6>(interface, fieldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeResidual<sphericalTensor8>(interface, fieldName, globalMaxJumpRes, globalMaxFluxRes);

    // maxTypeResidual<tensor2>(interface, fieldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeResidual<tensor4>(interface, fieldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeResidual<tensor6>(interface, fieldName, globalMaxJumpRes, globalMaxFluxRes);
    // maxTypeResidual<tensor8>(interface, fieldName, globalMaxJumpRes, globalMaxFluxRes);
}

template<class Type>
void Foam::dnaControl::writeResFields
(
    const regionInterface& interface,
    const word& fieldName,
    const Switch outputJumpResField,
    const Switch outputFluxResField,
    bool final
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;
    
    if (interface.meshA().foundObject<fieldType>(fieldName))
    {
        const fvPatchField<Type>& patchAField =
            interface.patchA().lookupPatchField<fieldType, Type>(fieldName);

        const fvPatchField<Type>& patchBField =
            interface.patchB().lookupPatchField<fieldType, Type>(fieldName);

        label patchAID = interface.patchAID();

        if 
        (
            isA<genericRegionCoupledJumpFvPatchField<Type>>(patchAField) &&
            isA<genericRegionCoupledFluxFvPatchField<Type>>(patchBField)
        )
        {
            const genericRegionCoupledJumpFvPatchField<Type>& jumpPatchAField = 
                refCast<const genericRegionCoupledJumpFvPatchField<Type>>(patchAField);

            const genericRegionCoupledFluxFvPatchField<Type>& fluxPatchBField = 
                refCast<const genericRegionCoupledFluxFvPatchField<Type>>(patchBField);

            if (outputJumpResField)
            {
                if (final)
                {
                    word jumpResFieldName = fieldName + "finalJumpRes";

                    finalJumpResFields_[jumpResFieldName]->boundaryField()[patchAID] = 
                        jumpPatchAField.rawResidual();
                }
                else
                {
                    word jumpResFieldName = fieldName + "initialJumpRes";

                    initialJumpResFields_[jumpResFieldName]->boundaryField()[patchAID] = 
                        jumpPatchAField.rawResidual();
                }
            }

            if (outputFluxResField)
            {
                if (final)
                {
                    word fluxResFieldName = fieldName + "finalFluxRes";

                    tmp<scalarField> resField =     
                        jumpPatchAField.interpolateFromNbrField
                        (
                            fluxPatchBField.rawResidual()
                        );
                    finalFluxResFields_[fluxResFieldName]->boundaryField()[patchAID] = resField();
                }
                else
                {
                    word fluxResFieldName = fieldName + "initialFluxRes";

                    tmp<scalarField> resField =     
                        jumpPatchAField.interpolateFromNbrField
                        (
                            fluxPatchBField.rawResidual()
                        );
                    initialFluxResFields_[fluxResFieldName]->boundaryField()[patchAID] = resField();
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
                refCast<const genericRegionCoupledFluxFvPatchField<Type>>(patchAField);

            const genericRegionCoupledJumpFvPatchField<Type>& jumpPatchBField = 
                refCast<const genericRegionCoupledJumpFvPatchField<Type>>(patchBField);

            if (outputFluxResField)
            {
                if (final)
                {
                    word fluxResFieldName = fieldName + "finalFluxRes";

                    finalFluxResFields_[fluxResFieldName]->boundaryField()[patchAID] = 
                        fluxPatchAField.rawResidual();
                }
                else
                {
                    word fluxResFieldName = fieldName + "initialFluxRes";

                    initialFluxResFields_[fluxResFieldName]->boundaryField()[patchAID] = 
                        fluxPatchAField.rawResidual();
                }
            }

            if (outputJumpResField)
            {
                if (final)
                {
                    word jumpResFieldName = fieldName + "finalJumpRes";

                    tmp<scalarField> resField =     
                        fluxPatchAField.interpolateFromNbrField
                        (
                            jumpPatchBField.rawResidual()
                        );
                    finalJumpResFields_[jumpResFieldName]->boundaryField()[patchAID] = resField();
                }
                else
                {
                    word jumpResFieldName = fieldName + "initialJumpRes";

                    tmp<scalarField> resField =     
                        fluxPatchAField.interpolateFromNbrField
                        (
                            jumpPatchBField.rawResidual()
                        );
                    initialJumpResFields_[jumpResFieldName]->boundaryField()[patchAID] = resField();
                }
            }
        }
    }
}

void Foam::dnaControl::writeResFields
(
    const regionInterface& interface,
    const word& fieldName,
    const Switch outputJumpResField,
    const Switch outputFluxResField,
    bool final
)
{
    writeResFields<scalar>(interface, fieldName, outputJumpResField, outputFluxResField, final);
    writeResFields<vector>(interface, fieldName, outputJumpResField, outputFluxResField, final);
    // writeResFields<sphericalTensor>(interface, fieldName, outputJumpResField, outputFluxResField, final);
    // writeResFields<symmTensor>(interface, fieldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFields<tensor>(interface, fieldName, outputJumpResFieldoutputFluxResField, final);

    // writeResFields<vector2>(interface, fieldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFields<vector4>(interface, fieldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFields<vector6>(interface, fieldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFields<vector8>(interface, fieldName, outputJumpResFieldoutputFluxResField, final);

    // writeResFields<sphericalTensor2>(interface, fieldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFields<sphericalTensor4>(interface, fieldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFields<sphericalTensor6>(interface, fieldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFields<sphericalTensor8>(interface, fieldName, outputJumpResFieldoutputFluxResField, final);

    // writeResFields<tensor2>(interface, fieldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFields<tensor4>(interface, fieldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFields<tensor6>(interface, fieldName, outputJumpResFieldoutputFluxResField, final);
    // writeResFields<tensor8>(interface, fieldName, outputJumpResFieldoutputFluxResField, final);
}

bool Foam::dnaControl::criteriaSatisfied()
{
    if (dnaResidualControl_.empty())
    {
        return false;
    }

    // no checks on first iteration - nothing has been calculated yet
    // only write initial residual fields
    if (corr_ == 1)
    {
        forAll (interfaces_, intI)
        {
            forAll(dnaResidualControl_, controlI)
            {
                word subFieldName = dnaResidualControl_[controlI].name;
                Switch outputJumpResField = dnaResidualControl_[controlI].outputJumpResField;
                Switch outputFluxResField = dnaResidualControl_[controlI].outputFluxResField;
                writeResFields(interfaces_[intI], subFieldName, outputJumpResField , outputFluxResField, false);
            }
        }
        return false;
    }

    //- List of max interface residuals for each subField
    List<scalar> globalMaxJumpRes(dnaResidualControl_.size());
    List<scalar> globalMaxFluxRes(dnaResidualControl_.size());
    globalMaxJumpRes = 0;
    globalMaxFluxRes = 0;

    forAll (interfaces_, intI)
    {
        forAll(dnaResidualControl_, controlI)
        {
            word subFieldName = dnaResidualControl_[controlI].name;
            maxResidual(interfaces_[intI], subFieldName, globalMaxJumpRes[controlI] , globalMaxFluxRes[controlI]);
        }
    }

    bool criteriaSatisfied = true;

    forAll(dnaResidualControl_, controlI)
    {
        criteriaSatisfied =
            criteriaSatisfied && (globalMaxJumpRes[controlI] <= dnaResidualControl_[controlI].maxJumpRes);
        
        criteriaSatisfied =
            criteriaSatisfied && (globalMaxFluxRes[controlI] <= dnaResidualControl_[controlI].maxFluxRes);
    }

    // output initial and final residual
    if(corr_ == 2)
    {
        forAll(dnaResidualControl_, controlI)
        {
            word fieldName = dnaResidualControl_[controlI].name;
            Info<< "    " << fieldName << " initialMaxJumpRes: " << globalMaxJumpRes[controlI] << nl 
                << "    " << fieldName << " initialMaxFluxRes: " << globalMaxFluxRes[controlI] << endl;
        }
    }

    if(criteriaSatisfied || corr_ == maxCoupleIter_)
    {
        forAll(dnaResidualControl_, controlI)
        {
            word fieldName = dnaResidualControl_[controlI].name;
            Info<< "    " << fieldName << " finalMaxJumpRes: " << globalMaxJumpRes[controlI] << nl 
                << "    " << fieldName << " finalMaxFluxRes: " << globalMaxFluxRes[controlI] << endl;
        }

        forAll (interfaces_, intI)
        {
            forAll(dnaResidualControl_, controlI)
            {
                word subFieldName = dnaResidualControl_[controlI].name;
                Switch outputJumpResField = dnaResidualControl_[controlI].outputJumpResField;
                Switch outputFluxResField = dnaResidualControl_[controlI].outputFluxResField;
                writeResFields(interfaces_[intI], subFieldName, outputJumpResField , outputFluxResField, true);
            }
        }
    }

    return criteriaSatisfied;

    // return false;
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dnaControl::createResFields()
{
    forAll (interfaces_, intI)
    {
        const regionInterface& interface = interfaces_[intI];

        forAll(dnaResidualControl_, controlI)
        {
            fieldData& control = dnaResidualControl_[controlI];
            
            if (control.outputJumpResField)
            {
                // initialJumpResFields
                if (interface.meshA().foundObject<volScalarField>(control.name + "initialJumpRes"))
                {
                    initialJumpResFields_.set
                    (
                        control.name + "initialJumpRes",
                        const_cast<volScalarField*>
                        (
                            &interface.meshA().lookupObject<volScalarField>(control.name + "initialJumpRes")
                        )
                    );
                }
                else
                {
                    Info << "DNA-Control: Creating field " << control.name << "initialJumpRes" << endl;
                    initialJumpResFields_.set
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
                if (interface.meshA().foundObject<volScalarField>(control.name + "finalJumpRes"))
                {
                    finalJumpResFields_.set
                    (
                        control.name + "finalJumpRes",
                        const_cast<volScalarField*>
                        (
                            &interface.meshA().lookupObject<volScalarField>(control.name + "finalJumpRes")
                        )
                    );
                }
                else
                {
                    Info << "DNA-Control: Creating field " << control.name << "finalJumpRes" << endl;
                    finalJumpResFields_.set
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
                if (interface.meshA().foundObject<volScalarField>(control.name + "initialFluxRes"))
                {
                    initialFluxResFields_.set
                    (
                        control.name + "initialFluxRes",
                        const_cast<volScalarField*>
                        (
                            &interface.meshA().lookupObject<volScalarField>(control.name + "initialFluxRes")
                        )
                    );
                }
                else
                {
                    Info << "DNA-Control: Creating field " << control.name << "initialFluxRes" << endl;
                    initialFluxResFields_.set
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
                if (interface.meshA().foundObject<volScalarField>(control.name + "finalFluxRes"))
                {
                    finalFluxResFields_.set
                    (
                        control.name + "finalFluxRes",
                        const_cast<volScalarField*>
                        (
                            &interface.meshA().lookupObject<volScalarField>(control.name + "finalFluxRes")
                        )
                    );
                }
                else
                {
                    Info << "DNA-Control: Creating field " << control.name << "finalFluxRes" << endl;
                    finalFluxResFields_.set
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

Foam::dnaControl::dnaControl(const Time& runTime, const word& fieldName, const regionInterfaceList& interfaces)
:
    IOobject
    (
        "dnaControl",
        runTime.timeName(),
        runTime
    ),
    fieldName_(fieldName),
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
    dnaResidualControl_(),
    initialJumpResFields_(),
    initialFluxResFields_(),
    finalJumpResFields_(),
    finalFluxResFields_(),
    maxCoupleIter_(0),
    corr_(0)
{
	read();

    createResFields();

	Info<< nl;
	if (dnaResidualControl_.empty())
	{
			Info<< "DNA-Control: no residual control data found. "
			<< "Calculations will employ " << maxCoupleIter_
			<< " outer coupling loops" << nl << endl;
	}
	else
	{
			Info<< "DNA-Control: max iterations = " << maxCoupleIter_
				<< endl;
			forAll(dnaResidualControl_, i)
			{
				Info<< "    field " << dnaResidualControl_[i].name << token::TAB
					<< ": maxJumpRes " << dnaResidualControl_[i].maxJumpRes
					<< ", maxFluxRes " << dnaResidualControl_[i].maxFluxRes
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

    if (criteriaSatisfied())
    {
        Info<< "DNA: converged in " << corr_ - 1
            << " iterations" << endl;

        corr_ = 0;
        return false;
    }

    if (corr_ > maxCoupleIter_)
    {
        if ((!dnaResidualControl_.empty()))
        {
            Info<< "DNA: not converged within "
                << maxCoupleIter_ << " iterations" << endl;
        }

        corr_ = 0;
        return false;
    }

    return true;
}
