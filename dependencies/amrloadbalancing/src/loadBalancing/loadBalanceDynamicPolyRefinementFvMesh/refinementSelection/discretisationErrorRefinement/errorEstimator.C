/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

Author
    Franjo Juretic, Creative Fields Ltd.  All rights reserved.
    Matthias Niethammer, TU Darmstadt.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "errorEstimator.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "linear.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

namespace errorEstimation
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

namespace
{

template <class Type, bool>
Type calculateExtrapolatedValueHelper
(
    const label cellI,
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const List<tmp<volVectorField>>& grad,
    const vector& deltaVec
)
{
    Type retVal;

    for(direction cmptI=0;cmptI<pTraits<Type>::nComponents;++cmptI)
    {
        retVal[cmptI] = field[cellI][cmptI] + (deltaVec & grad[cmptI]()[cellI]);
    }

    return retVal;
}

template <>
scalar calculateExtrapolatedValueHelper<scalar, true>
(
    const label cellI,
    const GeometricField<scalar, fvPatchField, volMesh>& field,
    const List<tmp<volVectorField>>& grad,
    const vector& deltaVec
)
{
    return (field[cellI] + (deltaVec & grad[0]()[cellI]));
}

template <class Type>
Type calculateExtrapolatedValue
(
    const label cellI,
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const List<tmp<volVectorField>>& grad,
    const vector& deltaVec
)
{
    return calculateExtrapolatedValueHelper
    <
        Type,
        pTraits<Type>::nComponents == 1
    >
    (
        cellI,
        field,
        grad,
        deltaVec
    );
}

}


template<class Type>
tmp<volScalarField> estimateDiscretisationError
(
    const GeometricField<Type, fvPatchField, volMesh>& field
)
{
    const auto faceValues = Foam::linearInterpolate(field);

    const auto& mesh = field.mesh();
    //const auto& cells = mesh.cells();

    //const vectorField& faceCentres = mesh.faceCentres();
    //const vectorField& cellCentres = mesh.cellCentres();

    tmp<volScalarField> errorFieldTmp
    (
        new volScalarField
        (
            IOobject
            (
                "discretisationError"+field.name(),
                field.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            field.dimensions()
        )
    );
    volScalarField& errorField = errorFieldTmp();
    errorField.internalField() = pTraits<scalar>::zero;

//    Info << "Calculating discretisation error for the field "
//        << field.name() << endl;

    List<tmp<volVectorField>> cmptGrad(pTraits<Type>::nComponents);
    for(direction cmptI=0;cmptI<pTraits<Type>::nComponents;++cmptI)
    {
        cmptGrad[cmptI] = fvc::grad(field.component(cmptI));
    }

    // Alternative face based implementation (MN 06/2021)
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const vectorField& C = mesh.C();
    const vectorField& Cf = mesh.Cf();

    forAll(faceValues(), faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        // face owner error
        const auto extrapolatedValueOwn =
                calculateExtrapolatedValue
                (
                    own,
                    field,
                    cmptGrad,
                    Cf[faceI]-C[own]
                );

        const scalar fErrorOwn =
            mag(faceValues()[faceI] - extrapolatedValueOwn);
        errorField[own] = Foam::max(errorField[own], fErrorOwn);

        // face neighbour error
        const auto extrapolatedValueNei =
                calculateExtrapolatedValue
                (
                    nei,
                    field,
                    cmptGrad,
                    Cf[faceI]-C[nei]
                );

        const scalar fErrorNei =
            mag(faceValues()[faceI] - extrapolatedValueNei);
        errorField[nei] = Foam::max(errorField[nei], fErrorNei);
    }

    forAll(faceValues().boundaryField(), patchI)
    {
        const auto& pfV = faceValues().boundaryField()[patchI];
        const labelUList& cellAddr = pfV.patch().faceCells();
        if (pfV.coupled())
        {
            forAll(pfV, faceI)
            {
                label own = cellAddr[faceI];
                const vectorField& Cf = pfV.patch().Cf();

                const auto extrapolatedValueOwn =
                    calculateExtrapolatedValue
                    (
                        own,
                        field,
                        cmptGrad,
                        Cf[faceI]-C[own]
                    );
                const scalar fErrorOwn =
                    mag(pfV[faceI] - extrapolatedValueOwn);
                errorField[own] = Foam::max(errorField[own], fErrorOwn);
            }
        }
        else
        {
            tmp<vectorField> pDeltas = pfV.patch().delta();
            forAll(pfV, faceI)
            {
                label own = cellAddr[faceI];

                const auto extrapolatedValueOwn =
                    calculateExtrapolatedValue
                    (
                        own,
                        field,
                        cmptGrad,
                        pDeltas()[faceI]
                    );

                const scalar fErrorOwn =
                    mag(pfV[faceI] - extrapolatedValueOwn);
                errorField[own] = Foam::max(errorField[own], fErrorOwn);
            }
        }
    }
//  End of alternative face based implementation

//    # ifdef USE_OMP
//    # pragma omp parallel for schedule(dynamic, 10)
//    # endif
//    forAll(cells, cellI)
//    {
//        scalar maxError = 0.0;

//        const point& cCentre = cellCentres[cellI];

//        for(const label& faceI : cells[cellI])
//        {
//            const vector& fCentre = faceCentres[faceI];

//            const auto extrapolatedValue =
//                calculateExtrapolatedValue
//                (
//                    cellI,
//                    field,
//                    cmptGrad,
//                    fCentre-cCentre
//                );

//            if( faceI < mesh.nInternalFaces() )
//            {
//                const scalar fError =
//                    mag(faceValues()[faceI] - extrapolatedValue);

//                maxError = Foam::max(maxError, fError);
//            }
//            else
//            {
//                const label patchID = mesh.boundaryMesh().whichPatch(faceI);
//                const auto patchRange = mesh.boundaryMesh().range(patchID);
//                const label fI = faceI - patchRange.first();

//                if( field.boundaryField()[patchID].size() > 0 )
//                {
//                    const scalar fError =
//                        mag
//                        (
//                            field.boundaryField()[patchID][fI] -
//                            extrapolatedValue
//                        );

//                    maxError = Foam::max(maxError, fError);
//                }
//            }
//        }

//        errorField[cellI] = maxError;
//    }

//    Info << "Finished estimating error for field " << field.name() << endl;

    return errorFieldTmp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace errorEstimation

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
