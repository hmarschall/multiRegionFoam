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

                data.append(fd);
            }
        }

        if (residualDict.found("maxJumpRes") && residualDict.found("maxFluxRes"))
        {
            fieldData fd;

            fd.name = fieldName_;
            fd.maxJumpRes = readScalar(residualDict.lookup("maxJumpRes"));
            fd.maxFluxRes = readScalar(residualDict.lookup("maxFluxRes"));

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

    // scalar globalMaxJumpRawRes = 0;
    // scalar globalMaxFluxRawRes = 0;

    // scalar globalMaxJumpNormRes = 0;
    // scalar globalMaxFluxNormRes = 0;
    
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
            //- Jump residual from patchA
            // scalar patchMaxJumpRawRes = 
            //     gMax(refCast<const genericRegionCoupledJumpFvPatchField<Type>>(patchAField).rawResidual());
            // globalMaxJumpRawRes = 
            //     max
            //     (
            //         globalMaxJumpRawRes,
            //         patchMaxJumpRawRes
            //     );
            
            // scalar patchMaxJumpNormRes = refCast<const genericRegionCoupledJumpFvPatchField<Type>>(patchAField).maxNormResidual();
            // globalMaxJumpNormRes = 
            //     max
            //     (
            //         globalMaxJumpNormRes,
            //         patchMaxJumpNormRes
            //     );
            
            scalar patchMaxJumpRes = refCast<const genericRegionCoupledJumpFvPatchField<Type>>(patchAField).ofNormResidual();
            globalMaxJumpRes = 
                max
                (
                    globalMaxJumpRes,
                    patchMaxJumpRes
                );
            
            //- Flux residual from patchB 
            // scalar patchMaxFluxRawRes = 
            //     gMax(refCast<const genericRegionCoupledFluxFvPatchField<Type>>(patchBField).rawResidual());
            // globalMaxFluxRawRes = 
            //     max
            //     (
            //         globalMaxFluxRawRes,
            //         patchMaxFluxRawRes
            //     );
            
            // scalar patchMaxFluxNormRes = refCast<const genericRegionCoupledFluxFvPatchField<Type>>(patchBField).maxNormResidual();
            // globalMaxFluxNormRes = 
            //     max
            //     (
            //         globalMaxFluxNormRes,
            //         patchMaxFluxNormRes
            //     );
            
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
            //- Jump residual from patchB
            // scalar patchMaxJumpRawRes = 
            //     gMax(refCast<const genericRegionCoupledJumpFvPatchField<Type>>(patchBField).rawResidual());
            // globalMaxJumpRawRes = 
            //     max
            //     (
            //         globalMaxJumpRawRes,
            //         patchMaxJumpRawRes
            //     );
            
            // scalar patchMaxJumpNormRes = refCast<const genericRegionCoupledJumpFvPatchField<Type>>(patchBField).maxNormResidual();
            // globalMaxJumpNormRes = 
            //     max
            //     (
            //         globalMaxJumpNormRes,
            //         patchMaxJumpNormRes
            //     );
            
            scalar patchMaxJumpRes = refCast<const genericRegionCoupledJumpFvPatchField<Type>>(patchBField).ofNormResidual();
            globalMaxJumpRes = 
                max
                (
                    globalMaxJumpRes,
                    patchMaxJumpRes
                );
            
            //- Flux residual from patchA 
            // scalar patchMaxFluxRawRes = 
            //     gMax(refCast<const genericRegionCoupledFluxFvPatchField<Type>>(patchAField).rawResidual());
            // globalMaxFluxRawRes = 
            //     max
            //     (
            //         globalMaxFluxRawRes,
            //         patchMaxFluxRawRes
            //     );
            
            // scalar patchMaxFluxNormRes = refCast<const genericRegionCoupledFluxFvPatchField<Type>>(patchAField).maxNormResidual();
            // globalMaxFluxNormRes = 
            //     max
            //     (
            //         globalMaxFluxNormRes,
            //         patchMaxFluxNormRes
            //     );
            
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
        // << "    maxJumpRawRes: " << globalMaxJumpRawRes << nl 
        // << "    maxFluxRawRes: " << globalMaxJumpRawRes << nl
        // << "    maxJumpNormRes: " << globalMaxJumpNormRes << nl 
        // << "    maxFluxNormRes: " << globalMaxFluxNormRes << nl
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

bool Foam::dnaControl::criteriaSatisfied()
{
    // no checks on first iteration - nothing has been calculated yet
    if (corr_ == 1 || dnaResidualControl_.empty())
    {
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
            Info<< "    initialMaxJumpRes: " << globalMaxJumpRes[controlI] << nl 
                << "    initialMaxFluxRes: " << globalMaxFluxRes[controlI] << endl;
        }
    }

    if(criteriaSatisfied || corr_ == maxCoupleIter_)
    {
        forAll(dnaResidualControl_, controlI)
        {
            Info<< "    finalMaxJumpRes: " << globalMaxJumpRes[controlI] << nl 
                << "    finalMaxFluxRes: " << globalMaxFluxRes[controlI] << endl;
        }
    }

    return criteriaSatisfied;

    // return false;
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
    maxCoupleIter_(0),
    corr_(0)
{
	read();

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
